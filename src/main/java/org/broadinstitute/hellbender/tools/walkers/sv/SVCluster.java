package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.*;

import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

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
 *         Matching SV type, with the exception of DEL/DUP combinations when --enable-cnv is used.
 *     </li>
 *     <li>
 *         Matching break-end strands.
 *     </li>
 *     <li>
 *         Interval reciprocal overlap (inapplicable for BNDs).
 *     </li>
 *     <li>
 *         Distance between event end-points (break-end window).
 *     </li>
 *     <li>
 *         Sample reciprocal overlap, based on carrier status determined by available genotypes (GT fields). If
 *         no GT fields are called for a given variant, the tool attempts to find carriers based on copy number (CN fields)
 *         and sample ploidy (as determined by the number of entries in GT). Note that carrier status can only be determined
 *         from CN for DEL and DUP records.
 *     </li>
 * </ul>
 *
 * <p>For CNV defragmentation (DEFRAGMENT_CNV DEFRAGMENT_CNV algorithm), the tool uses single-linkage clustering based
 * on the following criteria:</p>
 * <ul>
 *     <li>
 *         Must be a DEL/DUP/CNV and only be supported by a depth algorithm.
 *     </li>
 *     <li>
 *         Matching SV type
 *     </li>
 *     <li>
 *         Overlapping and meeting sample reciprocal overlap (see above) after padding both sides of each variant by the
 *         specified fraction of event length specified by --defrag-padding-fraction.
 * </ul>
 *
 * <p>Interval overlap, break-end window, and sample overlap parameters are defined for three combinations of event types
 * using the ALGORITHMS field:</p>
 * <ul>
 *     <li>
 *         Depth-only - both variants have only "depth" ALGORITHMS
 *     </li>
 *     <li>
 *         Evidenced/PESR - both variants at least one non-depth entry in ALGORITHMS
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
 *         SVLEN - variant length (-1 for BNDs and insertions of unknown length)
 *     </li>
 *     <li>
 *         STRANDS - break-end strands ("++", "+-", "-+", or "--")
 *     </li>
 *     <li>
 *         CHROM2 / END2 - end coordinates (BNDs only).
 *     </li>
 *     <li>
 *         ALGORITHMS - list of supporting tools/algorithms. These names may be arbitrary, except for depth-based callers,
 *         which should be listed as "depth".
 *     </li>
 * </ul>
 *
 * <p>The tool generates a new VCF with clusters collapsed into single representative records. By default, a MEMBERS field
 * is generated that lists the input variant IDs contained in that record's cluster.</p>
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
    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";
    public static final String ENABLE_CNV_LONG_NAME = "enable-cnv";
    public static final String DEFRAG_PADDING_FRACTION_LONG_NAME = "defrag-padding-fraction";
    public static final String CONVERT_INV_LONG_NAME = "convert-inv-to-bnd";
    public static final String ALGORITHM_LONG_NAME = "algorithm";
    public static final String FAST_MODE_LONG_NAME = "fast-mode";
    public static final String OMIT_MEMBERS_LONG_NAME = "omit-members";

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
    private String outputFile;

    @Argument(
            doc = "If supplied, generate variant IDs with this prefix",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    private String variantPrefix = null;

    /**
     * When enabled, DEL/DUP variants can be clustered together and the resulting record with have SVTYPE of CNV.
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
    private SVCollapser.BreakpointSummaryStrategy breakpointSummaryStrategy = SVCollapser.BreakpointSummaryStrategy.MEDIAN_START_MEDIAN_END;

    @Argument(fullName = DEFRAG_PADDING_FRACTION_LONG_NAME,
            doc = "Padding as a fraction of variant length (CNV defragmentation only)",
            optional = true
    )
    private double defragPaddingFraction = CNVDefragmenter.DEFAULT_PADDING_FRACTION;

    @Argument(fullName = ALGORITHM_LONG_NAME,
            doc = "Clustering algorithm",
            optional = true
    )
    private CLUSTER_ALGORITHM algorithm = CLUSTER_ALGORITHM.SINGLE_LINKAGE;

    @Argument(fullName = JointGermlineCNVSegmentation.DEFRAGMENTATION_PADDING_LONG_NAME,
            doc = "Extend events by this fraction on each side when determining overlap to merge",
            optional = true
    )
    private double defragmentationPadding = CNVDefragmenter.DEFAULT_PADDING_FRACTION;

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameters = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private SVClusterEngine<SVCallRecord> clusterEngine;
    private Set<String> samples;
    private String currentContig;
    private int numVariantsWritten = 0;

    private static final List<Allele> FAST_MODE_DEFAULT_NON_CALL_ALLELES = Arrays.asList(Allele.NO_CALL);

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        samples = getSamplesForVariants();

        if (algorithm == CLUSTER_ALGORITHM.DEFRAGMENT_CNV) {
            clusterEngine =  new CNVDefragmenter(dictionary, defragPaddingFraction, clusterParameters.getDepthParameters().getSampleOverlap());
        } else if (algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE || algorithm == CLUSTER_ALGORITHM.MAX_CLIQUE) {
            final SVCollapser collapser = new SVCollapser(breakpointSummaryStrategy);
            final LocatableClusterEngine.CLUSTERING_TYPE type = algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE ? LocatableClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE : LocatableClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
            clusterEngine = new SVClusterEngine<>(dictionary, type, enableCnv, collapser::collapse);
            clusterEngine.setDepthOnlyParams(clusterParameters.getDepthParameters());
            clusterEngine.setMixedParams(clusterParameters.getMixedParameters());
            clusterEngine.setEvidenceParams(clusterParameters.getPESRParameters());
        } else {
            throw new UnsupportedOperationException("Unsupported algorithm: " + algorithm.name());
        }

        writer = createVCFWriter(Paths.get(outputFile));
        writer.writeHeader(createHeader());
        currentContig = null;
    }

    @Override
    public Object onTraversalSuccess() {
        processClusters();
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
        SVCallRecord call = SVCallRecordUtils.create(variant);
        if (fastMode) {
            // Strip out no-call and ref genotypes to save memory and compute
            final GenotypesContext filteredGenotypes = GenotypesContext.copy(call.getGenotypes().stream()
                    .filter(SVCallRecordUtils::isAltGenotype)
                    .collect(Collectors.toList()));
            call = SVCallRecordUtils.copyCallWithNewGenotypes(call, filteredGenotypes);
        }

        // Flush clusters if we hit the next contig
        if (!call.getContigA().equals(currentContig)) {
            if (currentContig != null) {
                processClusters();
            }
            currentContig = call.getContigA();
        }

        // Add to clustering buffer
        if (convertInversions) {
            SVCallRecordUtils.convertInversionsToBreakends(call).forEachOrdered(clusterEngine::add);
        } else {
            clusterEngine.add(call);
        }
    }

    private void processClusters() {
        logger.info("Processing contig " + currentContig + "...");
        List<SVCallRecord> output = clusterEngine.getOutput();
        logger.info("Writing to file...");
        write(output);
        logger.info("Contig " + currentContig + " completed.");
    }

    private void write(final List<SVCallRecord> calls) {
        calls.stream().sorted(SVCallRecordUtils.getCallComparator(dictionary))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
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
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.END2_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Second position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, 1, VCFHeaderLineType.String, "Second contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.STRANDS_ATTRIBUTE, 1, VCFHeaderLineType.String, "First and second strands"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.ALGORITHMS_ATTRIBUTE, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Source algorithms"));
        if (!omitMembers) {
            header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Cluster variant ids"));
        }

        // Required format lines
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));

        return header;
    }

    /**
     * Build variant context variant context.
     *
     * @param call the call
     * @return the variant context
     */
    public VariantContext buildVariantContext(final SVCallRecord call) {
        // In case we're using fast mode
        final GenotypesContext filledGenotypes = SVCallRecordUtils.fillMissingSamplesWithGenotypes(call.getGenotypes(), samples, FAST_MODE_DEFAULT_NON_CALL_ALLELES, Collections.emptyMap());

        // Assign new variant ID
        final String newId = variantPrefix == null ? call.getId() : String.format("%s%08x", variantPrefix, numVariantsWritten++);

        // Build new variant
        final SVCallRecord finalCall = new SVCallRecord(newId, call.getContigA(), call.getPositionA(), call.getStrandA(), call.getContigB(),
                call.getPositionB(), call.getStrandB(), call.getType(), call.getLength(), call.getAlgorithms(), call.getAlleles(),
                filledGenotypes, call.getAttributes());
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(finalCall);
        if (omitMembers) {
            builder.rmAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        }
        return builder.make();
    }

}
