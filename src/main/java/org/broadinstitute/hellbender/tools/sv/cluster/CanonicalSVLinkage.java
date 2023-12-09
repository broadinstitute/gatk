package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.List;

/**
 * <p>Main class for SV clustering. Two items are clustered together if they:</p>
 * <ul>
 *   <li>Match event types</li>
 *   <li>Match strands</li>
 *   <li>Match contigs</li>
 *   <li>Meet minimum interval reciprocal overlap (unless inter-chromosomal or a BND).</li>
 *   <li>Meet minimum sample reciprocal overlap using carrier status. Genotypes (GT fields) are used to determine carrier
 *   status first. If not found and the event is of type DEL or DUP, then copy number (CN fields) are used instead.
 *   In the latter case, it is expected that ploidy can be determined from the number of entries in GT.</li>
 *   <li>Start and end coordinates are within a defined distance.</li>
 * </ul>
 *
 * <p>Interval overlap, sample overlap, and coordinate proximity parameters are defined separately for depth-only/depth-only,
 * depth-only/PESR, and PESR/PESR item pairs using the {@link ClusteringParameters} class. Note that all criteria must be met for
 * two candidate items to cluster, unless the comparison is depth-only/depth-only, in which case only one of
 * interval overlap or break-end proximity must be met. For insertions with undefined length (SVLEN less than 1), interval
 * overlap is tested assuming the given start position and a length of
 * {@value CanonicalSVLinkage#INSERTION_ASSUMED_LENGTH_FOR_OVERLAP}.</p>
 */
public class CanonicalSVLinkage<T extends SVCallRecord> extends SVClusterLinkage<T> {

    protected final SAMSequenceDictionary dictionary;
    protected final boolean clusterDelWithDup;  // DEL/DUP multi-allelic site clustering
    protected ClusteringParameters depthOnlyParams;
    protected ClusteringParameters mixedParams;
    protected ClusteringParameters evidenceParams;

    public static final int INSERTION_ASSUMED_LENGTH_FOR_OVERLAP = 50;
    public static final int INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY = 1;

    public static final double DEFAULT_RECIPROCAL_OVERLAP_DEPTH_ONLY = 0.8;
    public static final double DEFAULT_SIZE_SIMILARITY_DEPTH_ONLY = 0;
    public static final int DEFAULT_WINDOW_DEPTH_ONLY = 0;
    public static final double DEFAULT_SAMPLE_OVERLAP_DEPTH_ONLY = 0;

    public static final double DEFAULT_RECIPROCAL_OVERLAP_MIXED = 0.8;
    public static final double DEFAULT_SIZE_SIMILARITY_MIXED = 0;
    public static final int DEFAULT_WINDOW_MIXED = 1000;
    public static final double DEFAULT_SAMPLE_OVERLAP_MIXED = 0;

    public static final double DEFAULT_RECIPROCAL_OVERLAP_PESR = 0.5;
    public static final double DEFAULT_SIZE_SIMILARITY_PESR = 0;
    public static final int DEFAULT_WINDOW_PESR = 500;
    public static final double DEFAULT_SAMPLE_OVERLAP_PESR = 0;

    protected static final Logger logger = LogManager.getLogger(CanonicalSVLinkage.class);

    public static final ClusteringParameters DEFAULT_DEPTH_ONLY_PARAMS =
            ClusteringParameters.createDepthParameters(
                    DEFAULT_RECIPROCAL_OVERLAP_DEPTH_ONLY,
                    DEFAULT_SIZE_SIMILARITY_DEPTH_ONLY,
                    DEFAULT_WINDOW_DEPTH_ONLY,
                    DEFAULT_SAMPLE_OVERLAP_DEPTH_ONLY);
    public static final ClusteringParameters DEFAULT_MIXED_PARAMS =
            ClusteringParameters.createMixedParameters(
                    DEFAULT_RECIPROCAL_OVERLAP_MIXED,
                    DEFAULT_SIZE_SIMILARITY_MIXED,
                    DEFAULT_WINDOW_MIXED,
                    DEFAULT_SAMPLE_OVERLAP_MIXED);
    public static final ClusteringParameters DEFAULT_PESR_PARAMS =
            ClusteringParameters.createPesrParameters(
                    DEFAULT_RECIPROCAL_OVERLAP_PESR,
                    DEFAULT_SIZE_SIMILARITY_PESR,
                    DEFAULT_WINDOW_PESR,
                    DEFAULT_SAMPLE_OVERLAP_PESR);

    /**
     * Create a new engine
     * @param dictionary sequence dictionary pertaining to clustered records
     * @param clusterDelWithDup enables clustering of DEL/DUP variants into CNVs
     */
    public CanonicalSVLinkage(final SAMSequenceDictionary dictionary,
                              final boolean clusterDelWithDup) {
        Utils.nonNull(dictionary);
        this.dictionary = dictionary;
        this.depthOnlyParams = DEFAULT_DEPTH_ONLY_PARAMS;
        this.mixedParams = DEFAULT_MIXED_PARAMS;
        this.evidenceParams = DEFAULT_PESR_PARAMS;
        this.clusterDelWithDup = clusterDelWithDup;
    }

    @Override
    public boolean areClusterable(final SVCallRecord a, final SVCallRecord b) {
        if (!typesMatch(a, b)) {
            return false;
        }
        // Only require matching strands if BND or INV type
        if ((a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.BND
                || a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INV)
                && !strandsMatch(a, b)) {
            return false;
        }
        // Checks appropriate parameter set
        if (evidenceParams.isValidPair(a, b)) {
            return clusterTogetherWithParams(a, b, evidenceParams, dictionary);
        } else if (mixedParams.isValidPair(a, b)) {
            return clusterTogetherWithParams(a, b, mixedParams, dictionary);
        } else if (depthOnlyParams.isValidPair(a, b)) {
            return clusterTogetherWithParams(a, b, depthOnlyParams, dictionary);
        } else {
            return false;
        }
    }

    /**
     * Tests if SVTYPEs match, allowing for DEL/DUP/CNVs to match in some cases.
     */
    protected boolean typesMatch(final SVCallRecord a, final SVCallRecord b) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType aType = a.getType();
        final GATKSVVCFConstants.StructuralVariantAnnotationType bType = b.getType();
        if (aType == bType && a.getComplexSubtype() == b.getComplexSubtype()) {
            return true;
        }
        // Allow CNVs to cluster with both DELs and DUPs, but only allow DEL/DUP clustering if enabled
        if (a.isSimpleCNV() && b.isSimpleCNV()) {
            if (clusterDelWithDup || (aType == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV || bType == GATKSVVCFConstants.StructuralVariantAnnotationType.CNV)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Test if breakend strands match. True if one has null strands.
     */
    protected boolean strandsMatch(final SVCallRecord a, final SVCallRecord b) {
        if (a.nullStrands() || b.nullStrands()) {
            return true;
        }
        return a.getStrandA() == b.getStrandA() && a.getStrandB() == b.getStrandB();
    }

    /**
     * Helper function to test for clustering between two items given a particular set of parameters.
     * Note that the records have already been subjected to the checks in {@link #areClusterable(SVCallRecord, SVCallRecord)}.
     */
    private static boolean clusterTogetherWithParams(final SVCallRecord a, final SVCallRecord b,
                                                     final ClusteringParameters params,
                                                     final SAMSequenceDictionary dictionary) {
        // Contigs match
        if (!(a.getContigA().equals(b.getContigA()) && a.getContigB().equals(b.getContigB()))) {
            return false;
        }

        // If complex, test complex intervals
        if (a.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX
                && b.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX &&
                !testComplexIntervals(a, b, params.getReciprocalOverlap(), params.getSizeSimilarity(), params.getWindow(), dictionary)) {
            return false;
        }

        // Reciprocal overlap and size similarity
        // Check bypassed if both are inter-chromosomal
        final Boolean hasReciprocalOverlapAndSizeSimilarity;
        if (a.isIntrachromosomal()) {
            final boolean hasReciprocalOverlap = testReciprocalOverlap(a, b, params.getReciprocalOverlap());
            final boolean hasSizeSimilarity = testSizeSimilarity(a, b, params.getSizeSimilarity());
            hasReciprocalOverlapAndSizeSimilarity = hasReciprocalOverlap && hasSizeSimilarity;
            if (params.requiresOverlapAndProximity() && !hasReciprocalOverlapAndSizeSimilarity) {
                return false;
            }
        } else {
            hasReciprocalOverlapAndSizeSimilarity = null;
        }

        // Breakend proximity
        final boolean hasBreakendProximity = testBreakendProximity(a, b, params.getWindow(), dictionary);
        // Use short-circuiting statements since sample overlap is the least scalable / slowest check
        if (hasReciprocalOverlapAndSizeSimilarity == null) {
            return hasBreakendProximity && hasSampleOverlap(a, b, params.getSampleOverlap());
        } else {
            if (params.requiresOverlapAndProximity()) {
                return hasReciprocalOverlapAndSizeSimilarity && hasBreakendProximity && hasSampleOverlap(a, b, params.getSampleOverlap());
            } else {
                return (hasReciprocalOverlapAndSizeSimilarity || hasBreakendProximity) && hasSampleOverlap(a, b, params.getSampleOverlap());
            }
        }
    }

    /**
     * Performs overlap testing on each pair of complex intervals in two records, requiring each pair to be
     * sufficiently similar by reciprocal overlap, size similarity, and breakend proximity.
     */
    private static boolean testComplexIntervals(final SVCallRecord a, final SVCallRecord b, final double overlapThreshold,
                                                final double sizeSimilarityThreshold, final int window,
                                                final SAMSequenceDictionary dictionary) {
        final List<SVCallRecord.ComplexEventInterval> intervalsA = a.getComplexEventIntervals();
        final List<SVCallRecord.ComplexEventInterval> intervalsB = b.getComplexEventIntervals();
        if (intervalsA.size() != intervalsB.size()) {
            return false;
        }
        final Iterator<SVCallRecord.ComplexEventInterval> iterA = intervalsA.iterator();
        final Iterator<SVCallRecord.ComplexEventInterval> iterB = intervalsB.iterator();
        for (int i = 0; i < intervalsA.size(); i++) {
            final SVCallRecord.ComplexEventInterval cpxIntervalA = iterA.next();
            final SVCallRecord.ComplexEventInterval cpxIintervalB = iterB.next();
            if (cpxIntervalA.getIntervalType() != cpxIintervalB.getIntervalType()) {
                return false;
            }
            final SimpleInterval intervalA = cpxIntervalA.getInterval();
            final SimpleInterval intervalB = cpxIintervalB.getInterval();
            if (!(IntervalUtils.isReciprocalOverlap(intervalA, intervalB, overlapThreshold)
                    && testSizeSimilarity(intervalA.getLengthOnReference(), intervalB.getLengthOnReference(), sizeSimilarityThreshold)
                    && testBreakendProximity(new SimpleInterval(intervalA.getContig(), intervalA.getStart(), intervalA.getStart()),
                    new SimpleInterval(intervalA.getContig(), intervalA.getEnd(), intervalA.getEnd()),
                    new SimpleInterval(intervalB.getContig(), intervalB.getStart(), intervalB.getStart()),
                    new SimpleInterval(intervalB.getContig(), intervalB.getEnd(), intervalB.getEnd()), window, dictionary))) {
                return false;
            }
        }
        return true;
    }

    private static boolean testReciprocalOverlap(final SVCallRecord a, final SVCallRecord b, final double threshold) {
        final SimpleInterval intervalA = new SimpleInterval(a.getContigA(), a.getPositionA(), a.getPositionA() + getLength(a, INSERTION_ASSUMED_LENGTH_FOR_OVERLAP) - 1);
        final SimpleInterval intervalB = new SimpleInterval(b.getContigA(), b.getPositionA(), b.getPositionA() + getLength(b, INSERTION_ASSUMED_LENGTH_FOR_OVERLAP) - 1);
        return IntervalUtils.isReciprocalOverlap(intervalA, intervalB, threshold);
    }

    private static boolean testSizeSimilarity(final SVCallRecord a, final SVCallRecord b, final double threshold) {
        return testSizeSimilarity(getLength(a, INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY),
                getLength(b, INSERTION_ASSUMED_LENGTH_FOR_SIZE_SIMILARITY), threshold);
    }

    private static boolean testSizeSimilarity(final int lengthA, final int lengthB, final double threshold) {
        return Math.min(lengthA, lengthB) / (double) Math.max(lengthA, lengthB) >= threshold;
    }

    private static boolean testBreakendProximity(final SVCallRecord a, final SVCallRecord b, final int window,
                                                 final SAMSequenceDictionary dictionary) {
        return testBreakendProximity(a.getPositionAInterval(), a.getPositionBInterval(),
                b.getPositionAInterval(), b.getPositionBInterval(), window, dictionary);
    }

    private static boolean testBreakendProximity(final SimpleInterval intervalA1, final SimpleInterval intervalA2,
                                                 final SimpleInterval intervalB1, final SimpleInterval intervalB2,
                                                 final int window, final SAMSequenceDictionary dictionary) {
        final SimpleInterval intervalA1Padded = intervalA1.expandWithinContig(window, dictionary);
        final SimpleInterval intervalA2Padded = intervalA2.expandWithinContig(window, dictionary);
        if (intervalA1Padded == null) {
            logger.warn("Invalid start position " + intervalA1.getContig() + ":" + intervalA1.getStart() +
                    " - record will not be matched");
            return false;
        }
        if (intervalA2Padded == null) {
            logger.warn("Invalid end position " + intervalA2.getContig() + ":" + intervalA2.getStart() +
                    " - record will not be matched");
            return false;
        }
        return intervalA1Padded.overlaps(intervalB1) && intervalA2Padded.overlaps(intervalB2);
    }

    /**
     * Gets event length used for overlap testing.
     */
    private static int getLength(final SVCallRecord record, final int missingInsertionLength) {
        Utils.validate(record.isIntrachromosomal(), "Record must be intra-chromosomal");
        if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            return record.getLength() == null ? missingInsertionLength : Math.max(record.getLength(), 1);
        } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.BND) {
            return record.getPositionB() - record.getPositionA() + 1;
        } else if (record.getType() == GATKSVVCFConstants.StructuralVariantAnnotationType.CPX) {
            return getLengthComplex(record, missingInsertionLength);
        } else {
            // TODO lengths less than 1 shouldn't be valid
            return Math.max(record.getLength() == null ? 1 : record.getLength(), 1);
        }
    }

    /**
     * Returns length to be used for complex variant matching
     */
    private static int getLengthComplex(final SVCallRecord record, final int missingInsertionLength) {
        if (!record.isIntrachromosomal() || record.getPositionA() == record.getPositionB()) {
            // Insertion types use the sum of cpx intervals, or the "missing" length if there are none
            if (record.getComplexEventIntervals().isEmpty()) {
                return missingInsertionLength;
            } else {
                return record.getComplexEventIntervals().stream().mapToInt(SVCallRecord.ComplexEventInterval::getLengthOnReference).sum();
            }
        } else {
            // Intervaled types just use the difference in coordinates
            return record.getPositionB() - record.getPositionA() + 1;
        }
    }

    @Override
    public int getMaxClusterableStartingPosition(final SVCallRecord record) {
        return Math.max(
                getMaxClusterableStartingPositionWithParams(record, record.isDepthOnly() ? depthOnlyParams : evidenceParams, dictionary),
                getMaxClusterableStartingPositionWithParams(record, mixedParams, dictionary)
        );
    }

    /**
     * Returns max feasible start position of an item clusterable with the given record, given a set of clustering parameters.
     */
    private static int getMaxClusterableStartingPositionWithParams(final SVCallRecord call,
                                                                   final ClusteringParameters params,
                                                                   final SAMSequenceDictionary dictionary) {
        final String contig = call.getContigA();
        final int contigLength = dictionary.getSequence(contig).getSequenceLength();

        // Breakend proximity window
        final int maxPositionByWindow = Math.min(call.getPositionA() + params.getWindow(), contigLength);

        // Don't use overlap for inter-chromosomal events
        if (!call.isIntrachromosomal()) {
            return maxPositionByWindow;
        }

        // Reciprocal overlap window
        final int maxPositionByOverlap;
        final int maxPosition = (int) (call.getPositionA() + (1.0 - params.getReciprocalOverlap()) * getLength(call, INSERTION_ASSUMED_LENGTH_FOR_OVERLAP));
        maxPositionByOverlap = Math.min(maxPosition, contigLength);

        if (params.requiresOverlapAndProximity()) {
            return Math.min(maxPositionByOverlap, maxPositionByWindow);
        } else {
            return Math.max(maxPositionByOverlap, maxPositionByWindow);
        }
    }

    public final ClusteringParameters getDepthOnlyParams() {
        return depthOnlyParams;
    }
    public final ClusteringParameters getMixedParams() {
        return mixedParams;
    }
    public final ClusteringParameters getEvidenceParams() {
        return evidenceParams;
    }

    public final void setDepthOnlyParams(ClusteringParameters depthOnlyParams) { this.depthOnlyParams = depthOnlyParams; }

    public final void setMixedParams(ClusteringParameters mixedParams) {
        this.mixedParams = mixedParams;
    }

    public final void setEvidenceParams(ClusteringParameters evidenceParams) {
        this.evidenceParams = evidenceParams;
    }

}
