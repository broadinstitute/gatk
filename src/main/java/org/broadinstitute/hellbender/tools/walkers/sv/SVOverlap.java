package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.*;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.*;

/**
 * <p>Clusters structural variants based on coordinates, event type, and supporting algorithms. Primary use cases include:</p>
 * <ul>
 *     <li>
 *         Clustering SVs produced by multiple callers, based on interval overlap, breakpoint proximity, and sample overlap.
 *     </li>
 *     <li>
 *         Merging multiple SV VCFs with disjoint sets of samples and/or variants.
 *     </li>
 *     <li>
 *         Defragmentation of copy number variants produced with depth-based callers.
 *     </li>
 * </ul>
 *
 * <p>Clustering tasks can be accomplished using one of two algorithms. The SINGLE_LINKAGE algorithm produces clusters
 * for which all members cluster with <i>at least one</i> other member. The MAX_CLIQUE algorithm, however,
 * requires that all members cluster with <i>every other</i> member. The latter is in general non-polynomial in time and
 * space but implemented to minimize computations by traversing variants ordered by start position and efficiently
 * finalizing "active" clusters that are determined to be complete.</p>
 *
 * <p>The tool determines whether two given variants should cluster based following criteria:</p>
 * <ul>
 *     <li>
 *         Matching SV type. DEL and DUP are considered matching SV types if --enable-cnv is used and merged into
 *         a multi-allelic CNV type.
 *     </li>
 *     <li>
 *         Matching breakend strands (BND and INV only)
 *     </li>
 *     <li>
 *         Interval reciprocal overlap (inapplicable for BNDs).
 *     </li>
 *     <li>
 *         Distance between corresponding event breakends (breakend window).
 *     </li>
 *     <li>
 *         Sample reciprocal overlap, based on carrier status determined by available genotypes (GT fields). If
 *         no GT fields are called for a given variant, the tool attempts to find carriers based on copy number
 *         (CN field) and sample ploidy (as determined by the ECN FORMAT field).
 *     </li>
 * </ul>
 *
 * <p>For CNV defragmentation (DEFRAGMENT_CNV algorithm), the tool uses single-linkage clustering based
 * on the following linkage criteria:</p>
 * <ul>
 *     <li>
 *         Must be a DEL/DUP/CNV and only be supported by a depth algorithm.
 *     </li>
 *     <li>
 *         Matching SV type
 *     </li>
 *     <li>
 *         Overlapping after padding both sides of each variant by the fraction of event length specified by
 *         --defrag-padding-fraction.
 *     </li>
 *     <li>
 *         Sample overlap fraction (described above) specified by --defrag-sample-overlap.
 *     </li>
 * </ul>
 *
 * <p>Interval overlap, breakend window, and sample overlap parameters are defined for three combinations of event types
 * using the ALGORITHMS field, which describes the type of evidence that was used to call the variant:</p>
 * <ul>
 *     <li>
 *         Depth-only - both variants have solely "depth" ALGORITHMS
 *     </li>
 *     <li>
 *         PESR (paired-end/split-read) - both variants have at least one non-depth entry in ALGORITHMS
 *     </li>
 *     <li>
 *         Mixed - one variant is depth-only and the other is PESR
 *     </li>
 * </ul>
 *
 * <p>Users must supply one or more VCFs containing SVs with the following info fields:</p>
 *
 * <ul>
 *     <li>
 *         SVTYPE - event type (DEL, DUP, CNV, INS, INV, BND)
 *     </li>
 *     <li>
 *         SVLEN - variant length for INS only, if known
 *     </li>
 *     <li>
 *         STRANDS - breakend strands ("++", "+-", "-+", or "--") (BND and INV only)
 *     </li>
 *     <li>
 *         CHROM2 / END2 - mate coordinates (BND only)
 *     </li>
 *     <li>
 *         ALGORITHMS - list of supporting tools/algorithms. These names may be arbitrary, except for depth-based
 *         callers, which should be listed as "depth".
 *     </li>
 * </ul>
 *
 * <p>In addition, the following FORMAT fields must be defined:</p>
 *
 * <ul>
 *     <li>
 *         GT - genotype alleles
 *     </li>
 *     <li>
 *         ECN - expected copy number (i.e. ploidy)
 *     </li>
 *     <li>
 *         CN - genotype copy number (DEL/DUP/CNV only)
 *     </li>
 * </ul>
 *
 * <p>Note that for CNVs (DEL, DUP, multi-allelic CNV), GT alleles are set according to the CN/ECN fields. In
 * some cases, (e.g. diploid DUPs with CN 4), allele phasing cannot be determined unambiguously and GT is set with
 * no-call alleles.</p>
 *
 * <p>The tool generates a new VCF with clusters collapsed into single representative records. By default, a MEMBERS
 * field is generated that lists the input variant IDs contained in that record's cluster.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more SV VCFs
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Clustered VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCluster \
 *       -V variants.vcf.gz \
 *       -O clustered.vcf.gz \
 *       --algorithm SINGLE_LINKAGE
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Clusters structural variants",
        oneLineSummary = "Clusters structural variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVOverlap extends MultiVariantWalker {

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    @Argument(
            doc = "VCF containing SVs to annotate",
            fullName = "callset-vcf"
    )
    private GATKPath callsetVcf;


    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameterArgs = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private ReferenceSequenceFile reference;
    private VariantContextWriter writer;
    private Set<String> samples;
    private PartitionedSVClusterEngine<PartitionedSVCallRecord> engine;
    private PartitionedCallSetLinkage linkage;
    private Iterator<VariantContext> callsetVariants;
    private VariantContext currentCallsetVariant;
    private Comparator<VariantContext> variantComparator;
    private long nextItemId = 0L;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        dictionary = reference.getSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = getSamplesForVariants();

        linkage.setDepthOnlyParams(clusterParameterArgs.getDepthParameters());
        linkage.setMixedParams(clusterParameterArgs.getMixedParameters());
        linkage.setEvidenceParams(clusterParameterArgs.getPESRParameters());
        linkage = new PartitionedCallSetLinkage(dictionary);
        engine = new PartitionedSVClusterEngine<>(linkage, dictionary);

        variantComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        callsetVariants = new FeatureDataSource<VariantContext>(callsetVcf.toString()).iterator();

        writer = createVCFWriter(outputFile);
        writer.writeHeader(createHeader());
    }

    @Override
    public Object onTraversalSuccess() {
        flush(true);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        while (currentCallsetVariant == null && callsetVariants.hasNext()) {
            currentCallsetVariant = callsetVariants.next();
            if (variantComparator.compare(currentCallsetVariant, variant) <= 0) {
                engine.add(new PartitionedSVCallRecord(SVCallRecordUtils.create(currentCallsetVariant), true), nextItemId++);
                currentCallsetVariant = null;
            }
        }
        engine.add(new PartitionedSVCallRecord(SVCallRecordUtils.create(variant), false), nextItemId++);
        flush(false);
    }

    private void flush(final boolean force) {
        engine.flush(force).stream().forEach(this::processCluster);
    }

    private void processCluster(final PartitionedOutputCluster<PartitionedSVCallRecord> cluster) {
        // TODO : implement
    }

    private VCFHeader createHeader() {
        final VCFHeader header = getHeaderForVariants();
        header.addMetaDataLine(new VCFFormatHeaderLine("MGT", 1, VCFHeaderLineType.String, "Genotype from matched variant"));
        return header;
    }

    private static final class PartitionedCallSetLinkage extends CanonicalSVLinkage<PartitionedSVCallRecord> {

        public PartitionedCallSetLinkage(final SAMSequenceDictionary dictionary) {
            super(dictionary, false);
        }

        @Override
        public boolean areClusterable(final SVCallRecord a, final SVCallRecord b) {
            if (!algorithmsMatch(a, b)) {
                return false;
            }
            return super.areClusterable(a, b);
        }

        @Override
        protected boolean typesMatch(final SVCallRecord a, final SVCallRecord b) {
            final StructuralVariantType aType = a.getType();
            final StructuralVariantType bType = b.getType();
            if (aType == bType) {
                return true;
            }
            if (a.getType() == StructuralVariantType.BND) {
                return checkBnd(a, b);
            }
            if (b.getType() == StructuralVariantType.BND) {
                return checkBnd(b, a);
            }
            return super.typesMatch(a, b);
        }

        private boolean checkBnd(final SVCallRecord a, final SVCallRecord b) {
            if (b.isDepthOnly()) {
                return false;
            }
            return a.getStrandA() == b.getStrandA() && a.getStrandB() == b.getStrandB();
        }

        protected boolean algorithmsMatch(final SVCallRecord a, final SVCallRecord b) {
            final List<String> algorithmsA = a.getAlgorithms();
            final List<String> algorithmsB = b.getAlgorithms();
            if (algorithmsA.isEmpty() && algorithmsB.isEmpty()) {
                return true;
            }
            if (algorithmsA.size() == 1 &&  algorithmsB.size() == 1 && algorithmsA.get(0).equals(algorithmsB.get(0))) {
                return true;
            }
            final Set<String> setA = new HashSet<>(algorithmsA);
            if (setA.containsAll(algorithmsB)) {
                return true;
            }
            final Set<String> setB = new HashSet<>(algorithmsB);
            return setB.containsAll(setA);
        }
    }

}
