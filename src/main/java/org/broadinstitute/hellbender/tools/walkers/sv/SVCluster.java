package org.broadinstitute.hellbender.tools.walkers.sv;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.cluster.*;

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
public final class SVCluster extends SVClusterWalker {

    public static final String DEFRAG_PADDING_FRACTION_LONG_NAME = "defrag-padding-fraction";
    public static final String DEFRAG_SAMPLE_OVERLAP_LONG_NAME = "defrag-sample-overlap";

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

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameterArgs = new SVClusterEngineArgumentsCollection();

    protected SVClusterEngine clusterEngine;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        if (algorithm == CLUSTER_ALGORITHM.DEFRAGMENT_CNV) {
            clusterEngine = SVClusterEngineFactory.createCNVDefragmenter(dictionary, altAlleleSummaryStrategy,
                    reference, defragPaddingFraction, defragSampleOverlapFraction);
        } else if (algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE || algorithm == CLUSTER_ALGORITHM.MAX_CLIQUE) {
            final SVClusterEngine.CLUSTERING_TYPE type = algorithm == CLUSTER_ALGORITHM.SINGLE_LINKAGE ?
                    SVClusterEngine.CLUSTERING_TYPE.SINGLE_LINKAGE : SVClusterEngine.CLUSTERING_TYPE.MAX_CLIQUE;
            clusterEngine = SVClusterEngineFactory.createCanonical(type, breakpointSummaryStrategy,
                    altAlleleSummaryStrategy, dictionary, reference, enableCnv,
                    clusterParameterArgs.getDepthParameters(), clusterParameterArgs.getMixedParameters(),
                    clusterParameterArgs.getPESRParameters());
        } else {
            throw new IllegalArgumentException("Unsupported algorithm: " + algorithm.name());
        }
    }

    @Override
    public Object onTraversalSuccess() {
        clusterEngine.flush().stream().forEach(this::write);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
    }

    @Override
    public void applyRecord(final SVCallRecord record) {
        clusterEngine.addAndFlush(record).stream().forEach(this::write);
    }
}
