package org.broadinstitute.hellbender.tools.sv.cluster;

import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.function.BiPredicate;

/**
 * Stores clustering parameters for different combinations of supporting algorithm types (depth-only/depth-only,
 * depth-only/PESR, and PESR/PESR)
 */
public class ClusteringParameters {

    private final double reciprocalOverlap;  // minimum fractional reciprocal overlap of event intervals
    private final int window;  // maximum distance between variant end-points
    private final double sampleOverlap; // minimum fractional carrier sample overlap

    // if true, both reciprocal overlap and window criteria must be met
    // if false, reciprocal overlap and/or window criteria must be met
    private final boolean requiresOverlapAndProximity;

    // returns true if two given records are the correct type of pair for this parameter set
    private final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate;

    public ClusteringParameters(final double reciprocalOverlap, final int window, final double sampleOverlap,
                                final boolean overlapAndProximity, final BiPredicate<SVCallRecord, SVCallRecord> validRecordsPredicate) {
        this.reciprocalOverlap = reciprocalOverlap;
        this.window = window;
        this.sampleOverlap = sampleOverlap;
        this.requiresOverlapAndProximity = overlapAndProximity;
        this.validRecordsPredicate = validRecordsPredicate;
    }

    public double getReciprocalOverlap() {
        return reciprocalOverlap;
    }

    public int getWindow() {
        return window;
    }

    public double getSampleOverlap() {
        return sampleOverlap;
    }

    public boolean requiresOverlapAndProximity() {
        return requiresOverlapAndProximity;
    }

    public boolean isValidPair(final SVCallRecord a, final SVCallRecord b) {
        return validRecordsPredicate.test(a, b);
    }

    public static ClusteringParameters createDepthParameters(final double reciprocalOverlap, final int window, final double sampleOverlap) {
        return new ClusteringParameters(reciprocalOverlap, window, sampleOverlap, false, (a,b) -> a.isDepthOnly() && b.isDepthOnly());
    }

    public static ClusteringParameters createMixedParameters(final double reciprocalOverlap, final int window, final double sampleOverlap) {
        return new ClusteringParameters(reciprocalOverlap, window, sampleOverlap, true, (a,b) -> a.isDepthOnly() != b.isDepthOnly());
    }

    public static ClusteringParameters createPesrParameters(final double reciprocalOverlap, final int window, final double sampleOverlap) {
        return new ClusteringParameters(reciprocalOverlap, window, sampleOverlap, true, (a,b) -> !a.isDepthOnly() && !b.isDepthOnly());
    }
}
