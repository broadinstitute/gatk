package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME;
import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.FLAG_FIELD_LOGIC_LONG_NAME;

/***
 * Base class for tools that a simple interface for utilizing {@link SVClusterEngine}. It handles input/output easily,
 * including output sorting with spilling to disk to avoid excessive memory usage.
 */
public abstract class SVClusterWalker extends MultiVariantWalker {
    public static final String PLOIDY_TABLE_LONG_NAME = "ploidy-table";
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";
    public static final String ENABLE_CNV_LONG_NAME = "enable-cnv";
    public static final String ALGORITHM_LONG_NAME = "algorithm";
    public static final String FAST_MODE_LONG_NAME = "fast-mode";
    public static final String OMIT_MEMBERS_LONG_NAME = "omit-members";
    public static final String DEFAULT_NO_CALL_LONG_NAME = "default-no-call";
    public static final String MAX_RECORDS_IN_RAM_LONG_NAME = "max-records-in-ram";
    public static final String SITES_ONLY_LONG_NAME = "sites-only";
    public static final String LOW_MEM_LONG_NAME = "low-mem";

    /**
     * Genotype-object budget bounding the in-RAM window of genotype-heavy sort buffers in low-mem mode.
     * A record carries one genotype per sample, so at very large sample counts a flat
     * {@code --max-records-in-ram} window can exhaust the heap; bound by genotype VOLUME instead.
     */
    static final long LOW_MEM_GENOTYPE_BUFFER_BUDGET = 12_000_000L;
    /** Floor on the in-RAM record count so the number of spill files stays manageable. */
    static final int LOW_MEM_SORT_BUFFER_FLOOR = 50;

    /**
     * The enum Cluster algorithm.
     */
    public enum CLUSTER_ALGORITHM {
        /**
         * Defragment cnv cluster algorithm. Not supported with stratification.
         */
        DEFRAGMENT_CNV,
        /**
         * Single linkage cluster algorithm.
         */
        SINGLE_LINKAGE,
        /**
         * Max clique cluster algorithm.
         */
        MAX_CLIQUE
    }

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected GATKPath outputFile;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLE and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv). Note this is required unless using --" + SITES_ONLY_LONG_NAME + ".",
            fullName = PLOIDY_TABLE_LONG_NAME,
            mutex = {SITES_ONLY_LONG_NAME}
    )
    private GATKPath ploidyTablePath;

    @Argument(
            doc = "If supplied, generate variant IDs with this prefix",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    protected String variantPrefix = null;

    /**
     * When enabled, DEL and DUP variants will be clustered together. The resulting records with have an SVTYPE of CNV.
     */
    @Argument(
            doc = "Enable clustering DEL/DUP variants together as CNVs (does not apply to CNV defragmentation)",
            fullName = ENABLE_CNV_LONG_NAME,
            optional = true
    )
    protected boolean enableCnv = false;

    /**
     * Results in substantial space and time costs for large sample sets by clearing genotypes that are not needed for
     * clustering, but any associated annotation fields will be set to null in the output.
     */
    @Argument(
            doc = "Fast mode. Drops hom-ref and missing genotype fields and emits them as missing.",
            fullName = FAST_MODE_LONG_NAME,
            optional = true
    )
    protected boolean fastMode = false;

    @Argument(
            doc = "Omit cluster member ID annotations",
            fullName = OMIT_MEMBERS_LONG_NAME,
            optional = true
    )
    protected boolean omitMembers = false;

    @Argument(fullName = BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative value for a breakpoint cluster.",
            optional = true)
    protected CanonicalSVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy =
            CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE;

    @Argument(fullName = JointGermlineCNVSegmentation.ALT_ALLELE_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative alt allele for non-CNV biallelic sites with " +
                    "different subtypes.",
            optional = true)
    protected CanonicalSVCollapser.AltAlleleSummaryStrategy altAlleleSummaryStrategy =
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE;

    @Argument(fullName = FLAG_FIELD_LOGIC_LONG_NAME,
            doc = "Logic for collapsing Flag type INFO and FORMAT fields",
            optional = true)
    protected CanonicalSVCollapser.FlagFieldLogic flagFieldLogic = CanonicalSVCollapser.FlagFieldLogic.OR;

    @Argument(fullName = ALGORITHM_LONG_NAME,
            doc = "Clustering algorithm",
            optional = true
    )
    protected CLUSTER_ALGORITHM algorithm = CLUSTER_ALGORITHM.SINGLE_LINKAGE;

    @Argument(
            doc = "Drop all samples from input VCF(s) before clustering. The ploidy table may be omitted if using this " +
                    "option.",
            fullName = SITES_ONLY_LONG_NAME,
            mutex = {}
    )
    private boolean sitesOnly = false;

    /**
     * Default genotypes are assigned when they cannot be inferred from the inputs, such as when VCFs with different
     * variants and samples are provided.
     */
    @Argument(fullName = DEFAULT_NO_CALL_LONG_NAME,
            doc = "Default to no-call GT (e.g. ./.) instead of reference alleles (e.g. 0/0) when a genotype is not" +
                    " available",
            optional = true
    )
    protected boolean defaultNoCall = false;

    @Argument(fullName = MAX_RECORDS_IN_RAM_LONG_NAME,
            doc = "When writing VCF files that need to be sorted, this will specify the number of records stored in " +
            "RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a " +
            "VCF file, and increases the amount of RAM needed.",
            optional=true)
    public int maxRecordsInRam = 10000;

    @Argument(fullName = LOW_MEM_LONG_NAME,
            doc = "Enable low-memory two-pass mode. In pass one, genotypes are stripped from records before " +
                    "clustering to reduce memory usage. In pass two, the input VCFs are re-read to collect " +
                    "genotypes for each cluster. Produces byte-identical output to the default single-pass mode " +
                    "when sample overlap thresholds are at their defaults (0). Incompatible with non-zero " +
                    "sample overlap clustering parameters.",
            optional = true
    )
    protected boolean lowMem = false;

    protected SAMSequenceDictionary dictionary;
    protected ReferenceSequenceFile reference;
    protected PloidyTable ploidyTable;
    protected SortingCollection<VariantContext> sortingBuffer;
    protected VariantContextWriter writer;
    protected VCFHeader header;
    protected Set<String> samples;
    /** Count of pass-2 sites spilled to the completed-site buffer; for low-mem progress logging. */
    private long lowMemCompletedSites = 0;
    /**
     * Shared sample->index dictionary for compact pass-1 CNV copy-state ({@link CnvSampleCopyState}).
     * Non-null only in low-mem mode when CNV sample overlap is active (so CNV items can't be carrier-dropped
     * and would otherwise retain a genotype per sample). Built once in {@link #onTraversalStart}.
     */
    private CnvSampleCopyState.Dictionary cnvCopyStateDictionary;
    protected String currentContig;
    protected int numVariantsBuilt = 0;

    // ===== Low-memory two-pass mode shared state =====

    /** Collapser used for site-level collapse (pass 1) and genotype re-collapse finalization (pass 2). Set by subclass. */
    protected CanonicalSVCollapser lowMemCollapser;

    /**
     * Routing table: {@code routing.get(pass1RecordIdx)} is an int[] of PlannedSite sequence numbers
     * (indices into {@link #plannedSites}) that the record contributes to. Entries are null for records
     * that were not claimed by any cluster (should not occur in normal operation). Built during pass 1 by
     * {@link #lowMemCaptureCluster} and {@link #lowMemCaptureUnclustered}.
     *
     * <p>A record can map to multiple planned sites when using MAX_CLIQUE (OQ2: multi-membership is handled
     * correctly because routing is the authoritative source). The int[] encoding avoids O(N) String-keyed
     * HashMap entries that the old {@code syntheticIdToSiteIdx} required, reducing GC pressure.</p>
     *
     * <p>The ArrayList is nulled after pass 1 completes (just before pass 2 starts) so the routing
     * information (which is no longer needed after pass 2 setup) can be GC'd.</p>
     */
    private ArrayList<int[]> routing = new ArrayList<>();

    /**
     * Planned output sites in pass-1 write order. Entries are nulled out after the site is flushed so
     * the JVM can reclaim the whole PlannedSite shell (fixes the "3b lingering-shell" issue). Growing
     * as an ArrayList avoids any pre-size requirement (OQ1 resolved: size is not known in advance).
     */
    private final ArrayList<PlannedSite> plannedSites = new ArrayList<>();

    /** Sequential counter assigned to each record processed in pass 1. */
    private int pass1RecordIdx = 0;

    /**
     * Disk-backed, siteSeq-ordered buffer of completed sites' finalized records during pass 2. Created in
     * {@link #runLowMemFinalize} and drained (then cleaned up) there. Holds at most
     * {@code --max-records-in-ram} sites in memory and spills the rest to {@code tmpDir}, so the
     * completed-site memory is bounded regardless of head-of-line ordering.
     */
    private SortingCollection<SpilledSite> completedSiteBuffer;

    /**
     * Stripped pass-1 record that carries its own pass-1 sequential index. This replaces the old
     * {@code strippedRecordToPass1Index} {@link IdentityHashMap}, which pinned every stripped record
     * (and its retained genotypes) for the entire traversal. By stamping the index on the record
     * itself, the only references to a stripped record are inside the cluster engine's active window;
     * once a cluster is finalized the engine drops the record and it becomes eligible for GC, so peak
     * pass-1 genotype memory is bounded by the active clustering window rather than N.
     *
     * <p>The index is read back in {@link #lowMemCaptureCluster} / {@link #lowMemCaptureUnclustered}
     * when the engine emits the cluster. Reading rather than removing means MAX_CLIQUE multi-membership
     * (a record appearing in multiple emitted clusters) is handled for free.</p>
     *
     * <p>It also carries a precomputed {@code carrierCount} and overrides {@link #getCarrierCount()} so
     * the record can be fed to the engine with its genotype objects dropped entirely (empty genotypes)
     * when clustering will not touch them — see {@link #canUseGenotypeLightPass1Items()}. The collapser's
     * REPRESENTATIVE breakpoint tiebreaker reads the count via {@code getCarrierCount()}, so it remains
     * correct without the genotypes. This bounds the pass-1 active-window weight to O(window) instead of
     * O(window × carriers), which matters at common-variant hotspots with 250k samples.</p>
     */
    private static final class IndexedSVCallRecord extends SVCallRecord implements CnvCopyStateProvider {
        private final int pass1Index;
        private final int carrierCount;
        /** Compact CNV copy state, for genotype-light CNV items under active sample overlap; else null. */
        private final CnvSampleCopyState cnvCopyState;

        IndexedSVCallRecord(final SVCallRecord base, final GenotypesContext genotypes, final int pass1Index,
                            final int carrierCount, final CnvSampleCopyState cnvCopyState) {
            super(base.getId(), base.getContigA(), base.getPositionA(), base.getStrandA(), base.getContigB(),
                    base.getPositionB(), base.getStrandB(), base.getType(), base.getComplexSubtype(),
                    base.getComplexEventIntervals(), base.getLength(), base.getEvidence(), base.getAlgorithms(),
                    base.getAlleles(), genotypes, base.getAttributes(), base.getFilters(), base.getLog10PError());
            this.pass1Index = pass1Index;
            this.carrierCount = carrierCount;
            this.cnvCopyState = cnvCopyState;
        }

        int getPass1Index() {
            return pass1Index;
        }

        @Override
        public int getCarrierCount() {
            return carrierCount;
        }

        @Override
        public CnvSampleCopyState getCnvCopyState() {
            return cnvCopyState;
        }
    }

    /**
     * Container for a planned output site in low-memory two-pass mode. Holds the site-level record
     * (genotypes stripped) produced in pass 1 plus a running per-sample best-genotype accumulator
     * that is filled during pass 2. The entry in {@link #plannedSites} is nulled after flush so the
     * JVM can reclaim this entire shell (fixes the "3b lingering-shell" issue).
     */
    private static final class PlannedSite {
        /** Site-level record produced in pass 1 (genotypes list is empty / carrier-only). */
        final SVCallRecord siteRecord;
        /** Index into {@link #plannedSites} = pass-1 emission order = single-pass write order. */
        final int siteSeq;
        /**
         * If true, this site is an unclustered single-record passthrough (GroupedSVCluster's
         * no-stratum-match case). Its pass-2 genotypes are attached without collapse, matching
         * single-pass behaviour where the record is written as-is.
         */
        final boolean passthrough;
        /**
         * Running per-sample best genotype accumulated during pass 2, stored in a memory-compact form
         * ({@link FoldedGenotype}: serialized payload + cached comparator keys) rather than as a live
         * htsjdk {@link Genotype}, so the active-window footprint at large sample counts is bounded by
         * compact bytes instead of per-genotype attribute maps.
         */
        final Map<String, FoldedGenotype> bestGenotypes = new LinkedHashMap<>();
        /**
         * Number of pass-2 records still to be folded into this site. Set in pass 1 to the
         * member count; decremented as each member is read in pass 2. When it reaches 0 the
         * site is complete and can be finalized and freed.
         */
        int pendingMembers = 0;

        PlannedSite(final SVCallRecord siteRecord, final int siteSeq, final boolean passthrough) {
            this.siteRecord = siteRecord;
            this.siteSeq = siteSeq;
            this.passthrough = passthrough;
        }
    }

    /**
     * Memory-compact stand-in for an accumulated representative {@link Genotype} during pass-2 folding.
     *
     * <p>Holds the genotype serialized to a self-contained byte payload (all fields preserved, so the
     * finalized output genotype is byte-identical) plus the handful of scalar keys the representative-
     * genotype comparator chain reads. Storing bytes + keys instead of a live htsjdk {@link Genotype}
     * (object graph + extended-attribute map) shrinks the active-window footprint several-fold at large
     * sample counts, where one genotype per sample is held per open cluster.</p>
     *
     * <p>{@link #compare} replicates {@link CanonicalSVCollapser#getRepresentativeGenotype}'s comparator
     * chain exactly. The keys are extracted with the SAME {@code VariantContextGetters.getAttributeAsInt}
     * calls the comparators use, so the replica cannot drift from the source-of-truth chain. Any
     * divergence would surface in the differential unit test and the byte-identical integration tests.</p>
     */
    static final class FoldedGenotype {
        final byte[] payload;
        // Cached comparator keys (see the representative-genotype comparators in CanonicalSVCollapser).
        private final int nonRefCalledCount;  // # non-ref, called alleles
        private final int calledCount;        // # called alleles
        private final int quality;            // getAttributeAsInt(GQ, 0)
        private final int copyNumberQuality;  // getAttributeAsInt(CNQ, 0)
        private final int copyNumberDistance; // |ECN - CN(default 0)|
        private final boolean isDel;          // CN(default ECN) < ECN

        private FoldedGenotype(final byte[] payload, final int nonRefCalledCount, final int calledCount,
                               final int quality, final int copyNumberQuality,
                               final int copyNumberDistance, final boolean isDel) {
            this.payload = payload;
            this.nonRefCalledCount = nonRefCalledCount;
            this.calledCount = calledCount;
            this.quality = quality;
            this.copyNumberQuality = copyNumberQuality;
            this.copyNumberDistance = copyNumberDistance;
            this.isDel = isDel;
        }

        static FoldedGenotype of(final Genotype g) {
            int nonRefCalled = 0;
            int called = 0;
            for (final Allele a : g.getAlleles()) {
                if (a.isCalled()) {
                    called++;
                    if (a.isNonReference()) {
                        nonRefCalled++;
                    }
                }
            }
            final int quality = VariantContextGetters.getAttributeAsInt(g, VCFConstants.GENOTYPE_QUALITY_KEY, 0);
            final int cnq = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_QUALITY_FORMAT, 0);
            final int ecn = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, 0);
            final int cnForDistance = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, 0);
            final int cnForDel = VariantContextGetters.getAttributeAsInt(g, GATKSVVCFConstants.COPY_NUMBER_FORMAT, ecn);
            return new FoldedGenotype(encode(g), nonRefCalled, called, quality, cnq,
                    Math.abs(ecn - cnForDistance), cnForDel < ecn);
        }

        Genotype decode() {
            return decode(payload);
        }

        /**
         * Sign-identical to {@link CanonicalSVCollapser#getRepresentativeGenotype}'s comparator chain:
         * isNonRef, called, quality, nonRefCount (reversed: fewer wins), CN quality, CN distance
         * (reversed: closer to expected wins), DEL-over-DUP.
         */
        static int compare(final FoldedGenotype o1, final FoldedGenotype o2) {
            int c = Integer.compare(Math.min(1, o1.nonRefCalledCount), Math.min(1, o2.nonRefCalledCount));
            if (c != 0) { return c; }
            c = Integer.compare(o1.calledCount, o2.calledCount);
            if (c != 0) { return c; }
            c = Integer.compare(o1.quality, o2.quality);
            if (c != 0) { return c; }
            c = Integer.compare(o2.nonRefCalledCount, o1.nonRefCalledCount); // fewer non-ref alleles wins
            if (c != 0) { return c; }
            c = Integer.compare(o1.copyNumberQuality, o2.copyNumberQuality);
            if (c != 0) { return c; }
            c = Integer.compare(o2.copyNumberDistance, o1.copyNumberDistance); // closer to expected CN wins
            if (c != 0) { return c; }
            if (o1.isDel && !o2.isDel) { return 1; }
            if (o2.isDel && !o1.isDel) { return -1; }
            return 0;
        }

        // --- compact, self-contained genotype (de)serialization (DataOutput primitives, no object headers) ---

        private static byte[] encode(final Genotype g) {
            try {
                final ByteArrayOutputStream bos = new ByteArrayOutputStream(64);
                final DataOutputStream os = new DataOutputStream(bos);
                writeGenotypeTo(os, g);
                os.flush();
                return bos.toByteArray();
            } catch (final IOException e) {
                throw new UncheckedIOException("Failed to encode folded genotype", e);
            }
        }

        private static Genotype decode(final byte[] payload) {
            try {
                return readGenotypeFrom(new DataInputStream(new ByteArrayInputStream(payload)));
            } catch (final IOException e) {
                throw new UncheckedIOException("Failed to decode folded genotype", e);
            }
        }

        /**
         * Compact genotype serialization shared by {@link FoldedGenotype} and {@link SpilledSiteCodec}:
         * primitive {@link DataOutput} writes only (no {@link ObjectOutputStream}), so no per-stream object
         * handle table accumulates &mdash; the prior {@code ObjectOutputStream}-based site codec built a
         * multi-million-entry handle table when decoding an all-sample record, which is what OOM'd the drain.
         */
        static void writeGenotypeTo(final DataOutputStream os, final Genotype g) throws IOException {
            os.writeUTF(g.getSampleName());
            writeAlleles(os, g.getAlleles());
            os.writeBoolean(g.isPhased());
            os.writeInt(g.getGQ());
            os.writeInt(g.getDP());
            writeIntArray(os, g.hasAD() ? g.getAD() : null);
            writeIntArray(os, g.hasPL() ? g.getPL() : null);
            final String filters = g.getFilters();
            os.writeBoolean(filters != null);
            if (filters != null) {
                os.writeUTF(filters);
            }
            final Map<String, Object> ext = g.getExtendedAttributes();
            os.writeInt(ext.size());
            for (final Map.Entry<String, Object> e : ext.entrySet()) {
                os.writeUTF(e.getKey());
                writeAttributeValue(os, e.getValue());
            }
        }

        static Genotype readGenotypeFrom(final DataInputStream is) throws IOException {
            final String sample = is.readUTF();
            final List<Allele> alleles = readAlleles(is);
            final GenotypeBuilder gb = new GenotypeBuilder(sample, alleles);
            gb.phased(is.readBoolean());
            final int gq = is.readInt();
            if (gq >= 0) { gb.GQ(gq); }
            final int dp = is.readInt();
            if (dp >= 0) { gb.DP(dp); }
            final int[] ad = readIntArray(is);
            if (ad != null) { gb.AD(ad); }
            final int[] pl = readIntArray(is);
            if (pl != null) { gb.PL(pl); }
            if (is.readBoolean()) { gb.filters(is.readUTF()); }
            final int nAttr = is.readInt();
            final Map<String, Object> attrs = new LinkedHashMap<>();
            for (int i = 0; i < nAttr; i++) {
                final String key = is.readUTF();
                attrs.put(key, readAttributeValue(is));
            }
            gb.attributes(attrs);
            return gb.make();
        }

        static void writeAlleles(final DataOutputStream os, final List<Allele> alleles) throws IOException {
            os.writeInt(alleles.size());
            for (final Allele a : alleles) {
                os.writeBoolean(a.isNoCall());
                if (!a.isNoCall()) {
                    os.writeUTF(a.getDisplayString());
                    os.writeBoolean(a.isReference());
                }
            }
        }

        static List<Allele> readAlleles(final DataInputStream is) throws IOException {
            final int n = is.readInt();
            final List<Allele> alleles = new ArrayList<>(n);
            for (int i = 0; i < n; i++) {
                if (is.readBoolean()) {
                    alleles.add(Allele.NO_CALL);
                } else {
                    final String bases = is.readUTF();
                    alleles.add(Allele.create(bases, is.readBoolean()));
                }
            }
            return alleles;
        }

        private static void writeAttributeValue(final DataOutputStream os, final Object v) throws IOException {
            if (v == null) {
                os.writeByte(0);
            } else if (v instanceof Integer) {
                os.writeByte(1);
                os.writeInt((Integer) v);
            } else if (v instanceof String) {
                os.writeByte(2);
                os.writeUTF((String) v);
            } else if (v instanceof Double) {
                os.writeByte(3);
                os.writeDouble((Double) v);
            } else if (v instanceof Float) {
                os.writeByte(4);
                os.writeFloat((Float) v);
            } else if (v instanceof Long) {
                os.writeByte(5);
                os.writeLong((Long) v);
            } else if (v instanceof Boolean) {
                os.writeByte(6);
                os.writeBoolean((Boolean) v);
            } else if (v instanceof int[]) {
                os.writeByte(7);
                writeIntArray(os, (int[]) v);
            } else if (v instanceof List) {
                os.writeByte(8);
                final List<?> list = (List<?>) v;
                os.writeInt(list.size());
                for (final Object o : list) {
                    writeAttributeValue(os, o);
                }
            } else {
                // Fail loud rather than silently corrupt output: an unhandled FORMAT value type must be
                // added here explicitly (byte-identical output depends on preserving the runtime type).
                throw new IllegalStateException("Unsupported genotype FORMAT attribute type in low-mem fold: "
                        + v.getClass().getName());
            }
        }

        private static Object readAttributeValue(final DataInputStream is) throws IOException {
            final byte tag = is.readByte();
            switch (tag) {
                case 0: return null;
                case 1: return is.readInt();
                case 2: return is.readUTF();
                case 3: return is.readDouble();
                case 4: return is.readFloat();
                case 5: return is.readLong();
                case 6: return is.readBoolean();
                case 7: return readIntArray(is);
                case 8: {
                    final int n = is.readInt();
                    final List<Object> list = new ArrayList<>(n);
                    for (int i = 0; i < n; i++) {
                        list.add(readAttributeValue(is));
                    }
                    return list;
                }
                default:
                    throw new IllegalStateException("Unknown folded-genotype attribute tag: " + tag);
            }
        }

        private static void writeIntArray(final DataOutputStream os, final int[] arr) throws IOException {
            if (arr == null) {
                os.writeInt(-1);
            } else {
                os.writeInt(arr.length);
                for (final int v : arr) {
                    os.writeInt(v);
                }
            }
        }

        private static int[] readIntArray(final DataInputStream is) throws IOException {
            final int n = is.readInt();
            if (n < 0) {
                return null;
            }
            final int[] arr = new int[n];
            for (int i = 0; i < n; i++) {
                arr[i] = is.readInt();
            }
            return arr;
        }
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        dictionary = reference.getSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        // No REPRESENTATIVE-strategy restriction needed: pass-1 now retains carrier genotypes (non-CNV)
        // or CN/RD_CN/ECN-only genotypes (CNV), so the carrier-count tiebreaker used by REPRESENTATIVE
        // produces the same result as single-pass mode.
        if (ploidyTablePath == null) {
            ploidyTable = new PloidyTable(Collections.emptyMap());
        } else {
            ploidyTable = new PloidyTable(ploidyTablePath.toPath());
        }
        samples = sitesOnly ? Collections.emptySet() : getSamplesForVariants();
        // Shared sample->index dictionary for compact pass-1 CNV copy-state arrays (low-mem with active
        // CNV sample overlap only; otherwise unused). Built once so every compact item is index-aligned.
        cnvCopyStateDictionary = (lowMem && !canDropNonCarrierCnvGenotypes() && !samples.isEmpty())
                ? new CnvSampleCopyState.Dictionary(samples) : null;
        writer = createVCFWriter(outputFile);
        header = createHeader();
        writer.writeHeader(header);
        currentContig = null;
        final int sortBufferMaxInRam = lowMem ? genotypeHeavyMaxRecordsInRam() : maxRecordsInRam;
        if (lowMem) {
            logger.info(String.format(
                    "Low-mem output sort buffer: %d records in RAM (%d samples, --max-records-in-ram %d)",
                    sortBufferMaxInRam, samples.size(), maxRecordsInRam));
        }
        sortingBuffer = SortingCollection.newInstance(
                VariantContext.class,
                new VCFRecordCodec(header, true),
                header.getVCFRecordComparator(),
                sortBufferMaxInRam,
                tmpDir.toPath());
    }

    @Override
    public Object onTraversalSuccess() {
        for (final VariantContext variant : sortingBuffer) {
            writer.add(variant);
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (completedSiteBuffer != null) {
            completedSiteBuffer.cleanup();
        }
        if (sortingBuffer != null) {
            sortingBuffer.cleanup();
        }
        if (writer != null) {
            writer.close();
        }
    }

    /**
     * Subclasses should override this method
     */
    public abstract void applyRecord(final SVCallRecord record);

    /**
     * Whether it is safe to drop non-carrier genotypes from CNV records during low-mem pass 1. Only safe
     * when CNV sample-overlap linkage is disabled (all relevant sample-overlap thresholds are 0), since
     * {@code SVClusterLinkage.computeSampleOverlap} iterates every genotyped sample's copy state for CNV
     * records and would otherwise observe a different sample set. Defaults to false (retain all CNV
     * genotypes); subclasses override based on their clustering parameters.
     */
    protected boolean canDropNonCarrierCnvGenotypes() {
        return false;
    }

    /**
     * Whether low-mem pass 1 may feed the cluster engine genotype-free items (no genotype objects, only a
     * cached carrier count). Safe only when clustering never reads genotypes: the linkage is NOT
     * CNV-defragment (CNVLinkage reads carrier sets / copy state unconditionally) AND all sample-overlap
     * thresholds are 0 (so {@code computeSampleOverlap} is never invoked). Defaults to false (retain
     * genotypes); subclasses override based on their algorithm and clustering parameters.
     */
    protected boolean canUseGenotypeLightPass1Items() {
        return false;
    }

    // ===== Low-memory two-pass mode shared helpers =====

    /**
     * Pass 1 helper: strips genotypes from a record to the minimum needed for clustering, then
     * registers the stripped record for the two-pass fold. The original record ID is preserved so
     * MEMBERS annotations contain real input IDs. Object identity is tracked so the record can be
     * matched to its sequential pass-1 index during finalization.
     *
     * <p>Stripping strategy (mirrors fastMode to preserve sample-overlap clustering):
     * <ul>
     *   <li>Non-CNV: retain only carrier genotypes (same as {@code --fast-mode}).</li>
     *   <li>CNV: reduce every genotype to the CN/RD_CN/ECN attributes needed to determine carrier
     *       status (GT is zeroed). When {@link #canDropNonCarrierCnvGenotypes()} is true (CNV
     *       sample-overlap linkage disabled), non-carriers are additionally dropped — selected via
     *       {@link SVCallRecord#getCarrierGenotypeList()} on the reduced record, which uses the same
     *       isCarrier evaluation the engine/collapser apply, so the carrier set and count are
     *       preserved. Otherwise all reduced genotypes are kept, because
     *       {@code SVClusterLinkage.computeSampleOverlap} iterates EVERY genotyped sample's copy state
     *       for CNV records and would see a different sample set if non-carriers were removed.</li>
     * </ul>
     *
     * @param record full record delivered to {@code applyRecord}
     * @return genotype-stripped copy, carrying its pass-1 index, to feed to the cluster engine in pass 1
     */
    protected SVCallRecord lowMemStripAndRegister(final SVCallRecord record) {
        final GenotypesContext strippedGenotypes;
        final int carrierCount;
        CnvSampleCopyState cnvCopyState = null;
        if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            // CNV: reduce every genotype to CN/RD_CN/ECN with a no-call GT (the form the engine and
            // collapser will see). GT and all other FORMAT fields are dropped to save memory.
            final List<Genotype> cnOnly = record.getGenotypes().stream()
                    .map(g -> {
                        final GenotypeBuilder gb = new GenotypeBuilder(g.getSampleName());
                        // Preserve only the three copy-number attributes used by isCarrier / getCarrierSampleSet
                        final Object cn = g.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT);
                        final Object rdCn = g.getExtendedAttribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT);
                        final Object ecn = g.getExtendedAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT);
                        if (cn != null) gb.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, cn);
                        if (rdCn != null) gb.attribute(GATKSVVCFConstants.DEPTH_GENOTYPE_COPY_NUMBER_FORMAT, rdCn);
                        if (ecn != null) gb.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, ecn);
                        return gb.make();
                    })
                    .collect(Collectors.toList());
            final GenotypesContext cnReduced = GenotypesContext.create(new ArrayList<>(cnOnly));
            if (canDropNonCarrierCnvGenotypes()) {
                // CNV sample-overlap linkage is disabled, so computeSampleOverlap is never invoked.
                // Keep only carriers (selected on the reduced record so the carrier set / count match
                // exactly), dropping the non-carrier majority to bound pass-1 memory. The full CN-only
                // list is transient and freed once carriers are extracted.
                strippedGenotypes = GenotypesContext.copy(
                        SVCallRecordUtils.copyCallWithNewGenotypes(record, cnReduced).getCarrierGenotypeList());
                carrierCount = strippedGenotypes.size(); // all retained genotypes are carriers
            } else {
                // CNV sample-overlap may run; it iterates every genotyped sample's copy state. Rather than
                // RETAIN one Genotype per sample (~240k objects per item -> OOM in dense regions), capture
                // the per-sample copy state into compact int[] arrays and feed the engine a genotype-light
                // item. computeSampleOverlap reads the arrays via CnvCopyStateProvider, byte-identically.
                carrierCount = SVCallRecordUtils.copyCallWithNewGenotypes(record, cnReduced)
                        .getCarrierGenotypeList().size();
                cnvCopyState = CnvSampleCopyState.fromGenotypes(cnReduced, cnvCopyStateDictionary);
                strippedGenotypes = GenotypesContext.NO_GENOTYPES; // copy state now lives in cnvCopyState
            }
        } else {
            // Non-CNV: keep only carrier genotypes, exactly like --fast-mode
            strippedGenotypes = GenotypesContext.copy(record.getCarrierGenotypeList());
            carrierCount = strippedGenotypes.size(); // all retained genotypes are carriers
        }
        // When clustering will not read genotypes (non-defrag engine with sample-overlap disabled), feed
        // the engine an item with NO genotype objects to bound the pass-1 active window to O(window)
        // instead of O(window × carriers). The carrier count is preserved on the item for the
        // REPRESENTATIVE collapse tiebreaker (getCarrierCount), and linkage at sample-overlap 0 is
        // purely coordinate-based, so output is unchanged.
        final GenotypesContext engineGenotypes = canUseGenotypeLightPass1Items()
                ? GenotypesContext.NO_GENOTYPES : strippedGenotypes;
        // Stamp the pass-1 index on the record itself rather than pinning it in a global map, so the
        // record is GC-eligible as soon as the cluster engine finalizes its cluster.
        final SVCallRecord stripped =
                new IndexedSVCallRecord(record, engineGenotypes, pass1RecordIdx++, carrierCount, cnvCopyState);
        // Reserve a routing slot; actual site indices are filled in by lowMemCaptureCluster /
        // lowMemCaptureUnclustered when the cluster engine emits the cluster.
        routing.add(null);
        return stripped;
    }

    /**
     * Pass 1 collapser callback: collapses the cluster's site-level fields (no genotypes), appends
     * the site in write order, and maps each member's pass-1 index to the new site via the routing
     * table. The append order must match the order {@code write} would be called in single-pass mode
     * to keep output variant IDs byte-identical. Wire this as the {@code SVClusterEngine} collapser
     * function.
     *
     * <p>For MAX_CLIQUE (OQ2): a record can appear in multiple emitted clusters. Each appearance adds
     * the new site sequence number to the record's routing entry (int[]) via array extension. This is
     * safe because {@link #strippedRecordToPass1Index} is the lookup key and remains valid until
     * {@link #runLowMemFinalize} clears it before pass 2.</p>
     */
    protected SVCallRecord lowMemCaptureCluster(final SVClusterEngine.OutputCluster cluster) {
        final SVCallRecord siteRecord = lowMemCollapser.collapseWithGenotypes(cluster, Collections.emptyList());
        final int siteSeq = plannedSites.size();
        final PlannedSite site = new PlannedSite(siteRecord, siteSeq, false);
        plannedSites.add(site);
        int members = 0;
        for (final SVCallRecord item : cluster.getItems()) {
            final int p1Idx = ((IndexedSVCallRecord) item).getPass1Index();
            {
                // Register this site with the record's routing slot.
                final int[] existing = routing.get(p1Idx);
                if (existing == null) {
                    routing.set(p1Idx, new int[]{siteSeq});
                } else {
                    // MAX_CLIQUE multi-membership: append site to existing routing entry (OQ2 handled).
                    final int[] extended = Arrays.copyOf(existing, existing.length + 1);
                    extended[existing.length] = siteSeq;
                    routing.set(p1Idx, extended);
                }
                members++;
            }
        }
        site.pendingMembers = members;
        return siteRecord;
    }

    /**
     * Pass 1 helper for records that bypass clustering and are written as-is in single-pass mode
     * (e.g. {@code GroupedSVCluster} records that match no stratum). The supplied stripped record is
     * captured in write order as a passthrough site; its pass-2 genotypes are attached without any
     * collapse, matching single-pass behaviour. The {@code stripped} record must already have its
     * site-level annotations (MEMBERS, STRATUM) set by the caller.
     */
    protected void lowMemCaptureUnclustered(final SVCallRecord stripped) {
        final int siteSeq = plannedSites.size();
        final PlannedSite site = new PlannedSite(stripped, siteSeq, true);
        plannedSites.add(site);
        final int p1Idx = ((IndexedSVCallRecord) stripped).getPass1Index();
        routing.set(p1Idx, new int[]{siteSeq});
        site.pendingMembers = 1;
    }

    /**
     * Pass 2 + finalization for low-memory mode. Re-reads all input VCFs in the same coordinate-merged
     * order as pass 1 (via {@link MultiVariantDataSource}), folding each record's genotypes into the
     * accumulator of the site(s) it belongs to. Sites are written in pass-1 write order (so output IDs
     * match the single-pass baseline) via a write cursor that flushes and nulls each PlannedSite entry
     * as soon as it and all earlier sites are complete.
     *
     * <p>Memory-bounding behaviour: strippedRecordToPass1Index is cleared before pass 2 begins, freeing
     * all O(N) stripped records. The routing table is also freed (nulled) at that point. During pass 2
     * only the active write window is live in memory: sites with pendingMembers &gt; 0 retain their
     * bestGenotypes accumulator; completed sites are nulled in plannedSites[] immediately. Peak live
     * genotype memory is therefore bounded by the active clustering window, not by N or M.</p>
     *
     * <p>Call this from the subclass {@code onTraversalSuccess} after all cluster engine(s) have been
     * flushed.</p>
     */
    protected void runLowMemFinalize() {
        // Pass-1 stripped records are no longer reachable once the cluster engines have been flushed:
        // each carried its own pass-1 index (IndexedSVCallRecord) and was held only by the engine's
        // active window, so they have already become GC-eligible as their clusters were finalized.
        // The routing[] array (int[] per record) carries all remaining site-membership information.

        // Pass 2: fold genotypes
        // Re-read all input VCFs in the same coordinate-merged order as pass 1. Apply the walker's
        // traversal intervals so the pass-2 counter stays aligned with the pass-1 counter even when
        // the user supplied -L intervals (defect #3 fix).
        int pass2RecordIdx = 0;
        // Bounded, siteSeq-ordered buffer of completed sites' finalized records; spills to disk past
        // --max-records-in-ram so the completed-site memory is hard-bounded.
        completedSiteBuffer = SortingCollection.newInstance(
                SpilledSite.class,
                new SpilledSiteCodec(dictionary),
                Comparator.comparingInt(s -> s.siteSeq),
                maxRecordsInRam,
                tmpDir.toPath());
        try (final MultiVariantDataSource ds = new MultiVariantDataSource(getDrivingVariantsFeatureInputs(), 0)) {
            final List<SimpleInterval> traversalIntervals = getTraversalIntervals();
            if (traversalIntervals != null) {
                ds.setIntervalsForTraversal(traversalIntervals);
            }
            for (final VariantContext vc : ds) {
                final int p2Idx = pass2RecordIdx++;
                if (p2Idx > 0 && p2Idx % 200_000 == 0) {
                    final Runtime rt = Runtime.getRuntime();
                    logger.info(String.format(
                            "Low-mem pass 2: %d records folded, %d sites pending, %d completed; heap %d / %d MB",
                            p2Idx, plannedSites.size() - lowMemCompletedSites, lowMemCompletedSites,
                            (rt.totalMemory() - rt.freeMemory()) >> 20, rt.maxMemory() >> 20));
                }
                final int[] siteIndices = (p2Idx < routing.size()) ? routing.get(p2Idx) : null;
                if (siteIndices == null) {
                    continue;
                }
                // A record may map to multiple sites (MAX_CLIQUE multi-membership): fold into each.
                for (final int siteIdx : siteIndices) {
                    final PlannedSite site = plannedSites.get(siteIdx);
                    if (site == null) {
                        // Already flushed — this can happen if a site's last member was processed
                        // by a previous iteration and the site was flushed inline. Skip safely.
                        continue;
                    }
                    // Use the IDENTICAL per-record transform as apply() so that sitesOnly
                    // (zero genotypes) and fastMode (carrier-only) are both honoured in pass 2.
                    final SVCallRecord clusterRec = toClusterRecord(vc);
                    for (final Genotype g : clusterRec.getGenotypes()) {
                        final String sampleName = g.getSampleName();
                        final FoldedGenotype incoming = FoldedGenotype.of(g);
                        final FoldedGenotype current = site.bestGenotypes.get(sampleName);
                        // Mirror getRepresentativeGenotype(Arrays.asList(current, g)): Stream.max keeps the
                        // first maximal element, so the incoming genotype replaces the incumbent only when it
                        // strictly wins the comparator chain (ties keep the incumbent).
                        if (current == null || FoldedGenotype.compare(incoming, current) > 0) {
                            site.bestGenotypes.put(sampleName, incoming);
                        }
                    }
                    // Finalize the instant this site's last member is folded, spilling the whole
                    // record (site fields + genotypes) to the bounded completed-site buffer and freeing
                    // the plannedSites shell entirely. Spilling (rather than holding) caps live memory to
                    // the active clustering window plus the buffer's --max-records-in-ram window,
                    // independent of head-of-line order, and leaves no O(M) site records in RAM.
                    if (--site.pendingMembers <= 0) {
                        spillCompletedSite(site);
                        plannedSites.set(siteIdx, null); // release the whole shell (record + accumulator)
                    }
                }
            }
        }
        // Spill any sites whose members were not all re-read (edge cases in -L mode); normally none.
        for (int i = 0; i < plannedSites.size(); i++) {
            final PlannedSite ps = plannedSites.get(i);
            if (ps != null) {
                spillCompletedSite(ps);
                plannedSites.set(i, null);
            }
        }
        // Drain the completed-site buffer in siteSeq order (== single-pass write order, so output IDs
        // and same-coordinate tie-break ordering match the baseline) and write each finalized record.
        // The buffer holds only --max-records-in-ram sites in RAM; the rest are on disk.
        completedSiteBuffer.doneAdding();
        for (final SpilledSite spilled : completedSiteBuffer) {
            // The spilled record is the finalized site (site fields + genotypes), reconstructed by the
            // codec; write it directly. No plannedSites lookup is needed, so site records were freed at
            // completion time rather than held to O(M).
            write(SpilledSiteCodec.decodeRecord(spilled.payload, dictionary));
        }
        completedSiteBuffer.cleanup();
        completedSiteBuffer = null;
        // Release the routing table; it is no longer needed after pass 2 is complete.
        routing = null;
    }

    /**
     * In-RAM record count for genotype-heavy sort buffers in low-mem mode. Each output record carries
     * one genotype per sample, so a flat {@code --max-records-in-ram} window scales with sample count
     * and can exhaust the heap at large cohorts; bound the window by genotype VOLUME instead. The
     * external-sort OUTPUT is invariant to the spill count, so this stays byte-identical; at normal /
     * test sample counts the budget exceeds {@code --max-records-in-ram} and the count is unchanged.
     */
    private int genotypeHeavyMaxRecordsInRam() {
        return genotypeHeavyMaxRecordsInRam(maxRecordsInRam, samples == null ? 0 : samples.size());
    }

    /**
     * Pure arithmetic for {@link #genotypeHeavyMaxRecordsInRam()}, separated for unit testing. With no
     * samples, returns {@code maxRecordsInRam} unchanged; otherwise bounds the in-RAM record count so the
     * genotype VOLUME (records &times; samples) stays within {@link #LOW_MEM_GENOTYPE_BUFFER_BUDGET}, never
     * below {@link #LOW_MEM_SORT_BUFFER_FLOOR} and never above {@code maxRecordsInRam}.
     */
    static int genotypeHeavyMaxRecordsInRam(final int maxRecordsInRam, final int numSamples) {
        if (numSamples <= 0) {
            return maxRecordsInRam;
        }
        final long budgeted = LOW_MEM_GENOTYPE_BUFFER_BUDGET / numSamples;
        return (int) Math.min(maxRecordsInRam, Math.max(LOW_MEM_SORT_BUFFER_FLOOR, budgeted));
    }

    /**
     * Finalizes a completed planned site and appends the whole finalized record (site fields + genotypes)
     * to the bounded, siteSeq-ordered {@link #completedSiteBuffer} (which spills to disk past
     * {@code --max-records-in-ram}). The caller then nulls the {@link #plannedSites} entry, so neither the
     * accumulator nor the site record is retained in RAM after completion (the drain reconstructs the
     * record from the buffer).
     */
    private void spillCompletedSite(final PlannedSite site) {
        completedSiteBuffer.add(new SpilledSite(site.siteSeq,
                SpilledSiteCodec.encodeRecord(finalizePlannedSite(site))));
        lowMemCompletedSites++;
    }

    /**
     * Attaches the accumulated pass-2 genotypes to a planned site. Clustered sites run the same
     * post-collapse genotype processing as {@link CanonicalSVCollapser#collapse} (ref-allele
     * substitution then alt-allele harmonization); passthrough (unclustered) sites attach genotypes
     * raw, matching single-pass behaviour where the record is written without collapse.
     *
     * <p>This is the single shared helper for applying accumulated genotypes — used by both the
     * inline flush in {@link #runLowMemFinalize} and the safety-flush at the end (F1 stream parity:
     * one code path, no drift).</p>
     */
    private SVCallRecord finalizePlannedSite(final PlannedSite site) {
        final SVCallRecord rec = site.siteRecord;
        final List<Genotype> finalGenotypes;
        if (site.passthrough) {
            finalGenotypes = site.bestGenotypes.values().stream()
                    .map(FoldedGenotype::decode)
                    .collect(Collectors.toList());
        } else {
            final Allele refAllele = rec.getRefAllele();
            final List<Allele> altAlleles = rec.getAltAlleles();
            final List<Genotype> refSubstituted = site.bestGenotypes.values().stream()
                    .map(f -> lowMemCollapser.collapseSampleGenotypes(
                            Collections.singletonList(f.decode()), refAllele))
                    .collect(Collectors.toList());
            finalGenotypes = lowMemCollapser.harmonizeAltAlleles(altAlleles, refSubstituted);
        }
        return SVCallRecordUtils.copyCallWithNewGenotypes(rec, GenotypesContext.create(new ArrayList<>(finalGenotypes)));
    }

    /**
     * Shared per-record transform applied identically in pass 1 ({@link #apply}) and pass 2
     * ({@link #runLowMemFinalize}). Applies, in order:
     * <ol>
     *   <li>If {@code sitesOnly}: strip all genotypes so the resulting SVCallRecord has an empty
     *       genotypes list (no FORMAT column in output).</li>
     *   <li>Create the SVCallRecord from the (possibly genotype-stripped) VariantContext.</li>
     *   <li>If {@code fastMode} and the type is not CNV: reduce to carrier genotypes only.</li>
     * </ol>
     *
     * <p>By routing both passes through this single helper, sites-only and fast-mode transforms
     * can never drift between passes (F1 stream parity).</p>
     *
     * @param variant the raw VariantContext from the input VCF
     * @return the SVCallRecord that should be clustered (pass 1) or whose genotypes should be
     *         folded into the accumulator (pass 2)
     */
    private SVCallRecord toClusterRecord(VariantContext variant) {
        if (sitesOnly) {
            // Remove genotypes if we're in sites-only mode
            variant = new VariantContextBuilder(variant).genotypes(Collections.emptyList()).make();
        }
        SVCallRecord call = SVCallRecordUtils.create(variant, dictionary);
        if (fastMode && call.getType() != GATKSVVCFConstants.StructuralVariantAnnotationType.CNV) {
            // Strip out non-carrier genotypes to save memory and compute
            // Don't do for multi-allelic CNVs since carrier status can't be determined
            final GenotypesContext filteredGenotypes = GenotypesContext.copy(call.getCarrierGenotypeList());
            call = SVCallRecordUtils.copyCallWithNewGenotypes(call, filteredGenotypes);
        }
        return call;
    }

    @Override
    public void apply(VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final SVCallRecord call = toClusterRecord(variant);
        // Update current contig
        if (!call.getContigA().equals(currentContig)) {
            currentContig = call.getContigA();
            logger.info("Processing contig " + currentContig + "...");
        }
        applyRecord(call);
    }

    protected VCFHeader createHeader() {
        final VCFHeader header = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder(), samples);
        header.setSequenceDictionary(dictionary);

        // Required info lines
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVLEN));
        header.addMetaDataLine(GATKSVVCFHeaderLines.getInfoLine(GATKSVVCFConstants.SVTYPE));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END2_ATTRIBUTE, 1,
                VCFHeaderLineType.Integer, "Second position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, 1,
                VCFHeaderLineType.String, "Second contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.STRANDS_ATTRIBUTE, 1,
                VCFHeaderLineType.String, "First and second strands"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE,
                VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Source algorithms"));
        if (!omitMembers) {
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY,
                    VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Cluster variant ids"));
        }
        // Required format lines
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
        return header;
    }

    protected void write(final SVCallRecord call) {
        sortingBuffer.add(buildVariantContext(call));
    }

    protected VariantContext buildVariantContext(final SVCallRecord call) {
        // Add genotypes for missing samples
        final GenotypesContext filledGenotypes = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                call, samples, !defaultNoCall, ploidyTable, header);

        // Assign new variant ID
        final String newId = variantPrefix == null ? call.getId() : String.format("%s%08x", variantPrefix, numVariantsBuilt++);

        // Build new variant
        final SVCallRecord finalCall = new SVCallRecord(newId, call.getContigA(), call.getPositionA(), call.getStrandA(),
                call.getContigB(), call.getPositionB(), call.getStrandB(), call.getType(), call.getComplexSubtype(),
                call.getComplexEventIntervals(), call.getLength(), call.getEvidence(), call.getAlgorithms(), call.getAlleles(), filledGenotypes,
                call.getAttributes(), call.getFilters(), call.getLog10PError(), dictionary);
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(finalCall);
        if (omitMembers) {
            builder.rmAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        }
        return builder.make();
    }

    /**
     * A completed pass-2 site (finalized: site fields + genotypes), tagged with its siteSeq for ordered
     * draining. Package-private (not private) so {@code SpilledSiteCodec} round-trips can be unit-tested.
     */
    static final class SpilledSite {
        final int siteSeq;
        final byte[] payload;

        SpilledSite(final int siteSeq, final byte[] payload) {
            this.siteSeq = siteSeq;
            this.payload = payload;
        }
    }

    /**
     * Compact {@link SortingCollection} codec for {@link SpilledSite}. Serializes the whole finalized site
     * record using primitive {@link DataOutput} writes only &mdash; site-level fields inline, complex
     * intervals via their {@code encode()}/{@code decode()} string form, each genotype via the shared
     * {@link FoldedGenotype#writeGenotypeTo}/{@link FoldedGenotype#readGenotypeFrom} codec, and arbitrary
     * FORMAT/INFO value types via {@link FoldedGenotype#writeAttributeValue}. No {@link ObjectOutputStream}
     * is used, so decoding an all-sample record no longer builds a multi-million-entry object handle table
     * (the prior {@code ObjectInputStream}-based path OOM'd the drain at {@code HandleTable.grow}).
     * Genotypes are rebuilt with {@link GenotypeBuilder} so values round-trip exactly. Only the folded
     * (sparse) genotypes are stored, avoiding a full 250k-sample column per site. The sequence dictionary
     * is needed to reconstruct the record (coordinate validation + CPX intervals).
     *
     * <p>Package-private (not private) so the encode/decode round-trip can be unit-tested directly.</p>
     */
    static final class SpilledSiteCodec implements SortingCollection.Codec<SpilledSite> {
        private final SAMSequenceDictionary dictionary;
        private DataOutputStream out;
        private DataInputStream in;

        SpilledSiteCodec(final SAMSequenceDictionary dictionary) {
            this.dictionary = dictionary;
        }

        @Override
        public void setOutputStream(final OutputStream os) {
            out = new DataOutputStream(os);
        }

        @Override
        public void setInputStream(final InputStream is) {
            in = new DataInputStream(is);
        }

        @Override
        public void encode(final SpilledSite site) {
            try {
                // siteSeq + length-prefixed payload framing on the shared spill stream.
                out.writeInt(site.siteSeq);
                out.writeInt(site.payload.length);
                out.write(site.payload);
                out.flush();
            } catch (final IOException e) {
                throw new UncheckedIOException("Failed to encode spilled site", e);
            }
        }

        /**
         * Serializes a finalized site record to a self-contained byte payload. Done at add-time so the
         * completed-site buffer's in-RAM window holds compact bytes rather than fully-decoded genotype
         * graphs (each record carries up to one genotype per sample).
         */
        static byte[] encodeRecord(final SVCallRecord r) {
            try {
                final ByteArrayOutputStream bos = new ByteArrayOutputStream();
                final DataOutputStream os = new DataOutputStream(bos);
                os.writeUTF(r.getId());
                os.writeUTF(r.getContigA());
                os.writeInt(r.getPositionA());
                writeNullableBoolean(os, r.getStrandA());
                os.writeUTF(r.getContigB());
                os.writeInt(r.getPositionB());
                writeNullableBoolean(os, r.getStrandB());
                writeNullableString(os, r.getType() == null ? null : r.getType().name());
                writeNullableString(os, r.getComplexSubtype() == null ? null : r.getComplexSubtype().name());
                final List<SVCallRecord.ComplexEventInterval> cpx = r.getComplexEventIntervals();
                os.writeInt(cpx.size());
                for (final SVCallRecord.ComplexEventInterval ci : cpx) {
                    os.writeUTF(ci.encode());
                }
                writeNullableInt(os, r.getLength());
                final List<GATKSVVCFConstants.EvidenceTypes> evidence = r.getEvidence();
                os.writeInt(evidence.size());
                for (final GATKSVVCFConstants.EvidenceTypes e : evidence) {
                    os.writeUTF(e.name());
                }
                final List<String> algorithms = r.getAlgorithms();
                os.writeInt(algorithms.size());
                for (final String a : algorithms) {
                    os.writeUTF(a);
                }
                FoldedGenotype.writeAlleles(os, r.getAlleles());
                final Map<String, Object> attributes = r.getAttributes();
                os.writeInt(attributes.size());
                for (final Map.Entry<String, Object> e : attributes.entrySet()) {
                    os.writeUTF(e.getKey());
                    FoldedGenotype.writeAttributeValue(os, e.getValue());
                }
                final Set<String> filters = r.getFilters();
                os.writeInt(filters.size());
                for (final String f : filters) {
                    os.writeUTF(f);
                }
                writeNullableDouble(os, r.getLog10PError());
                final List<Genotype> genotypes = r.getGenotypes();
                os.writeInt(genotypes.size());
                for (final Genotype g : genotypes) {
                    FoldedGenotype.writeGenotypeTo(os, g);
                }
                os.flush();
                return bos.toByteArray();
            } catch (final IOException e) {
                throw new UncheckedIOException("Failed to encode spilled site", e);
            }
        }

        @Override
        public SpilledSite decode() {
            final int siteSeq;
            try {
                siteSeq = in.readInt();
            } catch (final EOFException e) {
                return null; // end of stream
            } catch (final IOException e) {
                throw new UncheckedIOException("Failed to decode spilled site", e);
            }
            try {
                final int len = in.readInt();
                final byte[] payload = new byte[len];
                in.readFully(payload);
                return new SpilledSite(siteSeq, payload);
            } catch (final IOException e) {
                throw new UncheckedIOException("Failed to decode spilled site", e);
            }
        }

        /** Reconstructs a finalized site record from a payload produced by {@link #encodeRecord}. */
        static SVCallRecord decodeRecord(final byte[] payload, final SAMSequenceDictionary dictionary) {
            try {
                final DataInputStream is = new DataInputStream(new ByteArrayInputStream(payload));
                final String id = is.readUTF();
                final String contigA = is.readUTF();
                final int positionA = is.readInt();
                final Boolean strandA = readNullableBoolean(is);
                final String contigB = is.readUTF();
                final int positionB = is.readInt();
                final Boolean strandB = readNullableBoolean(is);
                final String typeName = readNullableString(is);
                final GATKSVVCFConstants.StructuralVariantAnnotationType type = typeName == null ? null
                        : GATKSVVCFConstants.StructuralVariantAnnotationType.valueOf(typeName);
                final String cpxName = readNullableString(is);
                final GATKSVVCFConstants.ComplexVariantSubtype cpxSubtype = cpxName == null ? null
                        : GATKSVVCFConstants.ComplexVariantSubtype.valueOf(cpxName);
                final int nCpx = is.readInt();
                final List<SVCallRecord.ComplexEventInterval> cpx = new ArrayList<>(nCpx);
                for (int i = 0; i < nCpx; i++) {
                    cpx.add(SVCallRecord.ComplexEventInterval.decode(is.readUTF(), dictionary));
                }
                final Integer length = readNullableInt(is);
                final int nEvidence = is.readInt();
                final List<GATKSVVCFConstants.EvidenceTypes> evidence = new ArrayList<>(nEvidence);
                for (int i = 0; i < nEvidence; i++) {
                    evidence.add(GATKSVVCFConstants.EvidenceTypes.valueOf(is.readUTF()));
                }
                final int nAlgorithms = is.readInt();
                final List<String> algorithms = new ArrayList<>(nAlgorithms);
                for (int i = 0; i < nAlgorithms; i++) {
                    algorithms.add(is.readUTF());
                }
                final List<Allele> alleles = FoldedGenotype.readAlleles(is);
                final int nAttr = is.readInt();
                final Map<String, Object> attributes = new LinkedHashMap<>();
                for (int i = 0; i < nAttr; i++) {
                    final String key = is.readUTF();
                    attributes.put(key, FoldedGenotype.readAttributeValue(is));
                }
                final int nFilters = is.readInt();
                final Set<String> filters = new LinkedHashSet<>();
                for (int i = 0; i < nFilters; i++) {
                    filters.add(is.readUTF());
                }
                final Double log10PError = readNullableDouble(is);
                final int n = is.readInt();
                final List<Genotype> genotypes = new ArrayList<>(n);
                for (int i = 0; i < n; i++) {
                    genotypes.add(FoldedGenotype.readGenotypeFrom(is));
                }
                return new SVCallRecord(id, contigA, positionA, strandA, contigB, positionB,
                        strandB, type, cpxSubtype, cpx, length, evidence, algorithms, alleles, genotypes,
                        attributes, filters, log10PError, dictionary);
            } catch (final IOException e) {
                throw new UncheckedIOException("Failed to decode spilled site", e);
            }
        }

        @Override
        public SortingCollection.Codec<SpilledSite> clone() {
            return new SpilledSiteCodec(dictionary);
        }

        private static void writeNullableBoolean(final DataOutputStream os, final Boolean v) throws IOException {
            os.writeByte(v == null ? 0 : (v ? 1 : 2));
        }

        private static Boolean readNullableBoolean(final DataInputStream is) throws IOException {
            final byte b = is.readByte();
            return b == 0 ? null : b == 1;
        }

        private static void writeNullableInt(final DataOutputStream os, final Integer v) throws IOException {
            os.writeBoolean(v != null);
            if (v != null) {
                os.writeInt(v);
            }
        }

        private static Integer readNullableInt(final DataInputStream is) throws IOException {
            return is.readBoolean() ? is.readInt() : null;
        }

        private static void writeNullableDouble(final DataOutputStream os, final Double v) throws IOException {
            os.writeBoolean(v != null);
            if (v != null) {
                os.writeDouble(v);
            }
        }

        private static Double readNullableDouble(final DataInputStream is) throws IOException {
            return is.readBoolean() ? is.readDouble() : null;
        }

        private static void writeNullableString(final DataOutputStream os, final String v) throws IOException {
            os.writeBoolean(v != null);
            if (v != null) {
                os.writeUTF(v);
            }
        }

        private static String readNullableString(final DataInputStream is) throws IOException {
            return is.readBoolean() ? is.readUTF() : null;
        }
    }

}
