package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.*;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.*;

import static org.broadinstitute.hellbender.tools.walkers.sv.JointGermlineCNVSegmentation.BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME;

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
public final class SVCluster extends MultiVariantWalker {
    public static final String PLOIDY_TABLE_LONG_NAME = "ploidy-table";
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";
    public static final String ENABLE_CNV_LONG_NAME = "enable-cnv";
    public static final String DEFRAG_PADDING_FRACTION_LONG_NAME = "defrag-padding-fraction";
    public static final String DEFRAG_SAMPLE_OVERLAP_LONG_NAME = "defrag-sample-overlap";
    public static final String CONVERT_INV_LONG_NAME = "convert-inv-to-bnd";
    public static final String ALGORITHM_LONG_NAME = "algorithm";
    public static final String FAST_MODE_LONG_NAME = "fast-mode";
    public static final String OMIT_MEMBERS_LONG_NAME = "omit-members";
    public static final String INSERTION_LENGTH_SUMMARY_STRATEGY_LONG_NAME = "insertion-length-summary-strategy";
    public static final String DEFAULT_NO_CALL_LONG_NAME = "default-no-call";

    /**
     * The enum Cluster algorithm.
     */
    enum CLUSTER_ALGORITHM {
        /**
         * Defragment cnv cluster algorithm.
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
    private GATKPath outputFile;

    /**
     * Expected format is tab-delimited and contains a header with the first column SAMPLE and remaining columns
     * contig names. Each row corresponds to a sample, with the sample ID in the first column and contig ploidy
     * integers in their respective columns.
     */
    @Argument(
            doc = "Sample ploidy table (.tsv)",
            fullName = PLOIDY_TABLE_LONG_NAME
    )
    private GATKPath ploidyTablePath;

    @Argument(
            doc = "If supplied, generate variant IDs with this prefix",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    private String variantPrefix = null;

    /**
     * When enabled, DEL and DUP variants will be clustered together. The resulting records with have an SVTYPE of CNV.
     */
    @Argument(
            doc = "Enable clustering DEL/DUP variants together as CNVs (does not apply to CNV defragmentation)",
            fullName = ENABLE_CNV_LONG_NAME,
            optional = true
    )
    private boolean enableCnv = false;

    /**
     * When enabled, INV records will be converted to a pairs of BNDs prior to clustering.
     */
    @Argument(
            doc = "Convert inversions to BND records",
            fullName = CONVERT_INV_LONG_NAME,
            optional = true
    )
    private boolean convertInversions = false;

    /**
     * Results in substantial space and time costs for large sample sets by clearing genotypes that are not needed for
     * clustering, but any associated annotation fields will be set to null in the output.
     */
    @Argument(
            doc = "Fast mode. Drops hom-ref and no-call genotype fields and emits them as no-calls.",
            fullName = FAST_MODE_LONG_NAME,
            optional = true
    )
    private boolean fastMode = false;

    @Argument(
            doc = "Omit cluster member ID annotations",
            fullName = OMIT_MEMBERS_LONG_NAME,
            optional = true
    )
    private boolean omitMembers = false;

    @Argument(fullName = BREAKPOINT_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative value for a breakpoint cluster.",
            optional = true)
    private CanonicalSVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy =
            CanonicalSVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;

    @Argument(fullName = JointGermlineCNVSegmentation.ALT_ALLELE_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative alt allele for non-CNV biallelic sites with " +
                    "different subtypes.",
            optional = true)
    private CanonicalSVCollapser.AltAlleleSummaryStrategy altAlleleSummaryStrategy =
            CanonicalSVCollapser.AltAlleleSummaryStrategy.COMMON_SUBTYPE;

    @Argument(fullName = INSERTION_LENGTH_SUMMARY_STRATEGY_LONG_NAME,
            doc = "Strategy to use for choosing a representative value for insertion length when unknown.",
            optional = true)
    private CanonicalSVCollapser.InsertionLengthSummaryStrategy insertionLengthSummaryStrategy =
            CanonicalSVCollapser.InsertionLengthSummaryStrategy.MEDIAN;

    @Argument(fullName = DEFRAG_PADDING_FRACTION_LONG_NAME,
            doc = "Padding as a fraction of variant length for CNV defragmentation mode.",
            optional = true
    )
    private double defragPaddingFraction = CNVLinkage.DEFAULT_PADDING_FRACTION;

    @Argument(fullName = DEFRAG_SAMPLE_OVERLAP_LONG_NAME,
            doc = "Minimum sample overlap fraction. Use instead of --depth-sample-overlap in CNV defragmentation mode.",
            optional = true
    )
    private double defragSampleOverlapFraction = CNVLinkage.DEFAULT_SAMPLE_OVERLAP;

    @Argument(fullName = ALGORITHM_LONG_NAME,
            doc = "Clustering algorithm",
            optional = true
    )
    private CLUSTER_ALGORITHM algorithm = CLUSTER_ALGORITHM.SINGLE_LINKAGE;

    /**
     * Default genotypes are assigned when they cannot be inferred from the inputs, such as when VCFs with different
     * variants and samples are provided.
     */
    @Argument(fullName = DEFAULT_NO_CALL_LONG_NAME,
            doc = "Default to no-call GT (e.g. ./.) instead of reference alleles (e.g. 0/0) when a genotype is not" +
                    " available",
            optional = true
    )
    private boolean defaultNoCall = false;

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameterArgs = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private ReferenceSequenceFile reference;
    private PloidyTable ploidyTable;
    private VariantContextWriter writer;
    private SVClusterEngine<SVCallRecord> clusterEngine;
    private Set<String> samples;
    private String currentContig;
    private int numVariantsBuilt = 0;

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
        ploidyTable = new PloidyTable(ploidyTablePath.toPath());
        samples = getSamplesForVariants();

        if (algorithm == CLUSTER_ALGORITHM.DEFRAGMENT_CNV) {
            clusterEngine = SVClusterEngineFactory.createCNVDefragmenter(dictionary, altAlleleSummaryStrategy,
                    reference, defragPaddingFraction, defragSampleOverlapFraction);
        } else if (algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE || algorithm == CLUSTER_ALGORITHM.MAX_CLIQUE) {
            final SVClusterEngine.CLUSTERING_TYPE type = algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE ?
                    SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE : SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
            clusterEngine = SVClusterEngineFactory.createCanonical(type, breakpointSummaryStrategy,
                    altAlleleSummaryStrategy, insertionLengthSummaryStrategy, dictionary, reference, enableCnv,
                    clusterParameterArgs.getDepthParameters(), clusterParameterArgs.getMixedParameters(),
                    clusterParameterArgs.getPESRParameters());
        } else {
            throw new IllegalArgumentException("Unsupported algorithm: " + algorithm.name());
        }

        writer = createVCFWriter(outputFile);
        writer.writeHeader(createHeader());
        currentContig = null;
    }

    @Override
    public Object onTraversalSuccess() {
        write(true);
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
        final SVCallRecord call = SVCallRecordUtils.create(variant);
        final SVCallRecord filteredCall;
        if (fastMode) {
            // Strip out non-carrier genotypes to save memory and compute
            final GenotypesContext filteredGenotypes = GenotypesContext.copy(call.getCarrierGenotypeList());
            filteredCall = SVCallRecordUtils.copyCallWithNewGenotypes(call, filteredGenotypes);
        } else {
            filteredCall = call;
        }

        // Update current contig
        if (!filteredCall.getContigA().equals(currentContig)) {
            currentContig = filteredCall.getContigA();
            logger.info("Processing contig " + currentContig + "...");
        }

        // Add to clustering buffer
        if (convertInversions) {
            SVCallRecordUtils.convertInversionsToBreakends(filteredCall, dictionary).forEachOrdered(clusterEngine::add);
        } else {
            clusterEngine.add(filteredCall);
        }

        write(false);
    }

    private void write(final boolean force) {
        final List<SVCallRecord> records = force ? clusterEngine.forceFlush() : clusterEngine.flush();
        records.stream().map(this::buildVariantContext).forEachOrdered(writer::add);
    }

    private VCFHeader createHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samples);
        header.setVCFHeaderVersion(VCFHeaderVersion.VCF4_2);
        header.setSequenceDictionary(dictionary);

        // Copy from inputs
        getHeaderForVariants().getFormatHeaderLines().forEach(header::addMetaDataLine);
        getHeaderForVariants().getInfoHeaderLines().forEach(header::addMetaDataLine);

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

    public VariantContext buildVariantContext(final SVCallRecord call) {
        // Add genotypes for missing samples
        final GenotypesContext filledGenotypes = SVCallRecordUtils.populateGenotypesForMissingSamplesWithAlleles(
                call, samples, !defaultNoCall, ploidyTable);

        // Assign new variant ID
        final String newId = variantPrefix == null ? call.getId() : String.format("%s%08x", variantPrefix, numVariantsBuilt++);

        // Build new variant
        final SVCallRecord finalCall = new SVCallRecord(newId, call.getContigA(), call.getPositionA(), call.getStrandA(),
                call.getContigB(), call.getPositionB(), call.getStrandB(), call.getType(), call.getLength(),
                call.getAlgorithms(), call.getAlleles(), filledGenotypes, call.getAttributes(), dictionary);
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(finalCall);
        if (omitMembers) {
            builder.rmAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        }
        return builder.make();
    }

}
