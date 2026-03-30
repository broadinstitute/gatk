package org.broadinstitute.hellbender.tools.sv.cluster;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.List;
import java.util.Map;

/**
 * Linkage for matching complex SV records across different complex subtypes.
 * Requires that records have the same SV type but different complex subtypes,
 * and always uses complex interval testing for clustering.
 */
public class CrossSubtypeLinkage<T extends SVCallRecord> extends CanonicalSVLinkage<T> {

    public CrossSubtypeLinkage(final SAMSequenceDictionary dictionary) {
        super(dictionary, false);
    }

    @Override
    protected boolean typesMatch(final SVCallRecord a, final SVCallRecord b) {
        final GATKSVVCFConstants.StructuralVariantAnnotationType aType = a.getType();
        final GATKSVVCFConstants.StructuralVariantAnnotationType bType = b.getType();
        return aType == bType && a.getComplexSubtype() != b.getComplexSubtype();
    }

    @Override
    protected CanonicalLinkageResult clusterTogetherWithParams(final SVCallRecord a, final SVCallRecord b,
                                                               final ClusteringParameters params) {
        // Contigs match
        if (!(a.getContigA().equals(b.getContigA()) && a.getContigB().equals(b.getContigB()))) {
            return new CanonicalLinkageResult(false);
        }

        // Test overall reciprocal overlap, size similarity, and breakpoint distance
        final Integer breakpointDistance1 = getFirstBreakpointProximity(a, b);
        final Integer breakpointDistance2 = getSecondBreakpointProximity(a, b);
        final Double reciprocalOverlap = computeReciprocalOverlap(a, b);
        final Double sizeSimilarity = computeSizeSimilarity(a, b);

        if (!testBreakendProximity(breakpointDistance1, breakpointDistance2, params.getWindow())) {
            return new CanonicalLinkageResult(false);
        }
        if (!testReciprocalOverlap(reciprocalOverlap, params.getReciprocalOverlap())) {
            return new CanonicalLinkageResult(false);
        }
        if (!testSizeSimilarity(sizeSimilarity, params.getSizeSimilarity())) {
            return new CanonicalLinkageResult(false);
        }

        // dDUP/dDUP_iDEL vs dupINV/INVdup cross-subtype matching
        if (isDdupVsDupInv(a, b)) {
            return testDdupVsDupInv(a, b, params);
        }

        // Insertion-type cross-subtype matching (dDUP, dDUP_iDEL, INS_iDEL)
        if (isInsertionTypePair(a, b)) {
            return testInsertionTypePair(a, b, params);
        }

        // INV-anchored cross-subtype matching (dupINV/INVdup vs dupINVdup, delINV/INVdel vs delINVdel)
        if (isInvAnchoredPair(a, b)) {
            return testInvAnchoredPair(a, b, params);
        }

        // dDUP_iDEL vs dupINVdel/delINVdup cross-subtype matching
        if (isDdupIdelVsDupInvDelFamily(a, b)) {
            return testDdupInvWithDel(a, b, params);
        }

        // all other cases return false
        return new CanonicalLinkageResult(false);
    }

    /**
     * Tests whether a dDUP and dupINV/INVdup pair match by comparing their DUP intervals and checking
     * that the appropriate INV interval breakpoint is proximal to the dDUP's SINK_POS.
     * For dupINV, the 2nd breakpoint (end) of the INV is checked; for INVdup, the 1st breakpoint (start).
     *
     * Returns a {@link CanonicalLinkageResult} where:
     * - breakpointDistance1/2: distance between the sink position and the INV breakpoint (both set to the same value)
     * - reciprocalOverlap/sizeSimilarity: metrics for the DUP interval overlap
     * - interval arrays: DUP interval overlap details
     */
    private static CanonicalLinkageResult testDdupVsDupInv(final SVCallRecord a, final SVCallRecord b,
                                            final ClusteringParameters params) {
        final SVCallRecord ddup = isDdupFamily(a) ? a : b;
        final SVCallRecord other = isDdupFamily(a) ? b : a;

        // Verify the dDUP contains an INV interval in CPX_INTERVALS
        if (findIntervalByType(ddup, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) == null) {
            return new CanonicalLinkageResult(false);
        }

        // Find DUP intervals in both records
        final SVCallRecord.ComplexEventInterval ddupDup = findIntervalByType(ddup, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);
        final SVCallRecord.ComplexEventInterval otherDup = findIntervalByType(other, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);
        if (ddupDup == null || otherDup == null) {
            return new CanonicalLinkageResult(false);
        }

        // Test DUP interval match
        final CanonicalLinkageResult dupResult = testIntervalPair(ddupDup, otherDup, params);
        if (!dupResult.getResult()) {
            return dupResult;
        }

        // Test that the appropriate breakpoint of the INV interval is proximal to dDUP SINK_POS:
        // For dupINV, use the 2nd breakpoint (end); for INVdup, use the 1st breakpoint (start)
        final SVCallRecord.ComplexEventInterval otherInv = findIntervalByType(other, GATKSVVCFConstants.StructuralVariantAnnotationType.INV);
        if (otherInv == null || !hasSinkAttributes(ddup)) {
            return new CanonicalLinkageResult(false);
        }
        final String sinkChrom = getSinkChrom(ddup);
        if (!otherInv.getContig().equals(sinkChrom)) {
            return new CanonicalLinkageResult(false);
        }
        final int sinkPos = getSinkPos(ddup);
        final int invBreakpoint = other.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.INVdup
                ? otherInv.getStart()
                : otherInv.getEnd();
        final Integer sinkToInvDistance = Math.abs(invBreakpoint - sinkPos);
        if (sinkToInvDistance > params.getWindow()) {
            return new CanonicalLinkageResult(false);
        }

        return new CanonicalLinkageResult(true,
                dupResult.getReciprocalOverlap(), dupResult.getSizeSimilarity(),
                sinkToInvDistance, sinkToInvDistance,
                new Double[]{dupResult.getReciprocalOverlap()},
                new Double[]{dupResult.getSizeSimilarity()},
                new Integer[]{dupResult.getBreakpointDistance1()},
                new Integer[]{dupResult.getBreakpointDistance2()});
    }

    private static boolean isDdupVsDupInv(final SVCallRecord a, final SVCallRecord b) {
        return (isDdupFamily(a) && isDupInvOrInvDup(b)) || (isDupInvOrInvDup(a) && isDdupFamily(b));
    }

    private static boolean isInsertionTypePair(final SVCallRecord a, final SVCallRecord b) {
        return (isDdupFamily(a) || isInsIdel(a)) && (isInsIdel(b) || isDdupFamily(b));
    }

    private static boolean isInsIdel(final SVCallRecord record) {
        return record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.INS_iDEL;
    }

    /**
     * Tests whether a pair of insertion-type CPX variants (dDUP, dDUP_iDEL, INS_iDEL) match by checking
     * that they either both have or both lack an INV interval, that their SINK coordinates are proximal,
     * and that their source (DUP or INS) intervals overlap.
     *
     * Returns a {@link CanonicalLinkageResult} where:
     * - breakpointDistance1/2: SINK coordinate distances
     * - reciprocalOverlap/sizeSimilarity: source interval overlap metrics
     * - interval arrays: source interval overlap details
     */
    private static CanonicalLinkageResult testInsertionTypePair(final SVCallRecord a, final SVCallRecord b,
                                                             final ClusteringParameters params) {
        // Check that both have or both lack an INV interval
        final boolean aHasInv = findIntervalByType(a, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) != null;
        final boolean bHasInv = findIntervalByType(b, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) != null;
        if (aHasInv != bHasInv) {
            return new CanonicalLinkageResult(false);
        }

        // Check SINK coordinate proximity
        if (!hasSinkAttributes(a) || !hasSinkAttributes(b)) {
            return new CanonicalLinkageResult(false);
        }
        final String sinkChromA = getSinkChrom(a);
        final String sinkChromB = getSinkChrom(b);
        if (!sinkChromA.equals(sinkChromB)) {
            return new CanonicalLinkageResult(false);
        }
        final Integer sinkStartDist = Math.abs(getSinkPos(a) - getSinkPos(b));
        final Integer sinkEndDist = Math.abs(getSinkEnd(a) - getSinkEnd(b));
        if (!testBreakendProximity(sinkStartDist, sinkEndDist, params.getWindow())) {
            return new CanonicalLinkageResult(false);
        }

        // Check DUP interval overlap (for INS_iDEL, use the INS interval instead of DUP)
        final SVCallRecord.ComplexEventInterval intervalA = isInsIdel(a)
                ? findIntervalByType(a, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)
                : findIntervalByType(a, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);
        final SVCallRecord.ComplexEventInterval intervalB = isInsIdel(b)
                ? findIntervalByType(b, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)
                : findIntervalByType(b, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);
        if (intervalA == null || intervalB == null) {
            return new CanonicalLinkageResult(false);
        }
        final CanonicalLinkageResult sourceResult = testIntervalPair(intervalA, intervalB, params);
        if (!sourceResult.getResult()) {
            return sourceResult;
        }

        // Initialize interval metric arrays with source (DUP/INS) interval results
        Double[] intervalOverlaps = new Double[]{sourceResult.getReciprocalOverlap()};
        Double[] intervalSizeSims = new Double[]{sourceResult.getSizeSimilarity()};
        Integer[] intervalBpDist1s = new Integer[]{sourceResult.getBreakpointDistance1()};
        Integer[] intervalBpDist2s = new Integer[]{sourceResult.getBreakpointDistance2()};

        // If both records have DEL intervals, check their overlap and overwrite arrays
        final SVCallRecord.ComplexEventInterval delA = findIntervalByType(a, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final SVCallRecord.ComplexEventInterval delB = findIntervalByType(b, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        if (delA != null && delB != null) {
            final CanonicalLinkageResult delResult = testIntervalPair(delA, delB, params);
            if (!delResult.getResult()) {
                return delResult;
            }
            intervalOverlaps = new Double[]{sourceResult.getReciprocalOverlap(), delResult.getReciprocalOverlap()};
            intervalSizeSims = new Double[]{sourceResult.getSizeSimilarity(), delResult.getSizeSimilarity()};
            intervalBpDist1s = new Integer[]{sourceResult.getBreakpointDistance1(), delResult.getBreakpointDistance1()};
            intervalBpDist2s = new Integer[]{sourceResult.getBreakpointDistance2(), delResult.getBreakpointDistance2()};
        }

        return new CanonicalLinkageResult(true,
                sourceResult.getReciprocalOverlap(), sourceResult.getSizeSimilarity(),
                sinkStartDist, sinkEndDist,
                intervalOverlaps, intervalSizeSims, intervalBpDist1s, intervalBpDist2s);
    }

    /**
     * Checks whether one record has 2 intervals and the other has 3, and both contain an INV.
     */
    private static boolean isInvAnchoredPair(final SVCallRecord a, final SVCallRecord b) {
        final int aSize = a.getComplexEventIntervals().size();
        final int bSize = b.getComplexEventIntervals().size();
        if (!((aSize == 2 && bSize == 3) || (aSize == 3 && bSize == 2))) {
            return false;
        }
        return findIntervalByType(a, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) != null
                && findIntervalByType(b, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) != null;
    }

    /**
     * Tests whether two INV-containing complex subtypes match by anchoring on the INV interval
     * and comparing the non-INV interval from the shorter (2-interval) record with the positionally
     * corresponding interval from the longer (3-interval) record.
     * If INV is at index 0: compare shorter[1] with longer[2]
     * If INV is at index 1: compare shorter[0] with longer[0]
     * ie. dupINV vs. dupINVdup, INVdel vs. delINVdel, INVdup vs. delINVdup, etc.
     */
    private static CanonicalLinkageResult testInvAnchoredPair(final SVCallRecord a, final SVCallRecord b,
                                                              final ClusteringParameters params) {
        final List<SVCallRecord.ComplexEventInterval> aIntervals = a.getComplexEventIntervals();
        final List<SVCallRecord.ComplexEventInterval> bIntervals = b.getComplexEventIntervals();
        final List<SVCallRecord.ComplexEventInterval> shorter = aIntervals.size() == 2 ? aIntervals : bIntervals;
        final List<SVCallRecord.ComplexEventInterval> longer = aIntervals.size() == 3 ? aIntervals : bIntervals;

        // Find INV position in the shorter record
        final int invIdx = shorter.get(0).getIntervalSVType() == GATKSVVCFConstants.StructuralVariantAnnotationType.INV ? 0 : 1;

        // Compare INV intervals (shorter's INV vs longer's middle element)
        final CanonicalLinkageResult invResult = testIntervalPair(shorter.get(invIdx), longer.get(1), params);
        if (!invResult.getResult()) {
            return invResult;
        }

        // Compare non-INV interval: INV at 0 -> shorter[1] vs longer[2]; INV at 1 -> shorter[0] vs longer[0]
        final SVCallRecord.ComplexEventInterval shorterNonInv = shorter.get(1 - invIdx);
        final SVCallRecord.ComplexEventInterval longerNonInv = longer.get(invIdx == 0 ? 2 : 0);
        if (shorterNonInv.getIntervalSVType() != longerNonInv.getIntervalSVType()) {
            return new CanonicalLinkageResult(false);
        }
        final CanonicalLinkageResult nonInvResult = testIntervalPair(shorterNonInv, longerNonInv, params);
        if (!nonInvResult.getResult()) {
            return nonInvResult;
        }

        return new CanonicalLinkageResult(true,
                invResult.getReciprocalOverlap(), invResult.getSizeSimilarity(),
                invResult.getBreakpointDistance1(), invResult.getBreakpointDistance2(),
                new Double[]{invResult.getReciprocalOverlap(), nonInvResult.getReciprocalOverlap()},
                new Double[]{invResult.getSizeSimilarity(), nonInvResult.getSizeSimilarity()},
                new Integer[]{invResult.getBreakpointDistance1(), nonInvResult.getBreakpointDistance1()},
                new Integer[]{invResult.getBreakpointDistance2(), nonInvResult.getBreakpointDistance2()});
    }

    /**
     * Tests whether two complex event intervals match by breakpoint proximity, reciprocal overlap,
     * and size similarity. Returns a passing result with metrics, or a failing result.
     */
    private static CanonicalLinkageResult testIntervalPair(final SVCallRecord.ComplexEventInterval a,
                                                           final SVCallRecord.ComplexEventInterval b,
                                                           final ClusteringParameters params) {
        final Integer bpDist1 = getFirstBreakpointProximity(a, b);
        final Integer bpDist2 = getSecondBreakpointProximity(a, b);
        if (!testBreakendProximity(bpDist1, bpDist2, params.getWindow())) {
            return new CanonicalLinkageResult(false);
        }
        final Double overlap = computeReciprocalOverlap(a.getInterval(), b.getInterval());
        if (!testReciprocalOverlap(overlap, params.getReciprocalOverlap())) {
            return new CanonicalLinkageResult(false);
        }
        final Double sizeSim = computeSizeSimilarity(a.getInterval().size(), b.getInterval().size());
        if (!testSizeSimilarity(sizeSim, params.getSizeSimilarity())) {
            return new CanonicalLinkageResult(false);
        }
        return new CanonicalLinkageResult(true, overlap, sizeSim, bpDist1, bpDist2);
    }

    private static boolean isDdupFamily(final SVCallRecord record) {
        return record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.dDUP
                || record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL;
    }

    private static boolean isDupInvOrInvDup(final SVCallRecord record) {
        return record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.dupINV
                || record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.INVdup;
    }

    private static boolean isDdupIdel(final SVCallRecord record) {
        return record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.dDUP_iDEL;
    }

    private static boolean isDupInvDelFamily(final SVCallRecord record) {
        return record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.dupINVdel
                || record.getComplexSubtype() == GATKSVVCFConstants.ComplexVariantSubtype.delINVdup;
    }

    private static boolean isDdupIdelVsDupInvDelFamily(final SVCallRecord a, final SVCallRecord b) {
        return (isDdupIdel(a) && isDupInvDelFamily(b)) || (isDupInvDelFamily(a) && isDdupIdel(b));
    }

    /**
     * Tests whether a dDUP_iDEL and a dupINVdel/delINVdup pair match by comparing their DEL and DUP
     * intervals and verifying that the dDUP_iDEL contains an INV interval.
     */
    private static CanonicalLinkageResult testDdupInvWithDel(final SVCallRecord a, final SVCallRecord b,
                                                             final ClusteringParameters params) {
        final SVCallRecord ddupIdel = isDdupIdel(a) ? a : b;
        final SVCallRecord other = isDdupIdel(a) ? b : a;

        // Verify the dDUP_iDEL contains an INV interval
        if (findIntervalByType(ddupIdel, GATKSVVCFConstants.StructuralVariantAnnotationType.INV) == null) {
            return new CanonicalLinkageResult(false);
        }

        // Compare DEL intervals
        final SVCallRecord.ComplexEventInterval delA = findIntervalByType(ddupIdel, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        final SVCallRecord.ComplexEventInterval delB = findIntervalByType(other, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        if (delA == null || delB == null) {
            return new CanonicalLinkageResult(false);
        }
        final CanonicalLinkageResult delResult = testIntervalPair(delA, delB, params);
        if (!delResult.getResult()) {
            return delResult;
        }

        // Compare DUP intervals
        final SVCallRecord.ComplexEventInterval dupA = findIntervalByType(ddupIdel, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);
        final SVCallRecord.ComplexEventInterval dupB = findIntervalByType(other, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP);
        if (dupA == null || dupB == null) {
            return new CanonicalLinkageResult(false);
        }
        final CanonicalLinkageResult dupResult = testIntervalPair(dupA, dupB, params);
        if (!dupResult.getResult()) {
            return dupResult;
        }

        return new CanonicalLinkageResult(true,
                dupResult.getReciprocalOverlap(), dupResult.getSizeSimilarity(),
                dupResult.getBreakpointDistance1(), dupResult.getBreakpointDistance2(),
                new Double[]{dupResult.getReciprocalOverlap(), delResult.getReciprocalOverlap()},
                new Double[]{dupResult.getSizeSimilarity(), delResult.getSizeSimilarity()},
                new Integer[]{dupResult.getBreakpointDistance1(), delResult.getBreakpointDistance1()},
                new Integer[]{dupResult.getBreakpointDistance2(), delResult.getBreakpointDistance2()});
    }

    private static SVCallRecord.ComplexEventInterval findIntervalByType(
            final SVCallRecord record, final GATKSVVCFConstants.StructuralVariantAnnotationType type) {
        final List<SVCallRecord.ComplexEventInterval> intervals = record.getComplexEventIntervals();
        for (final SVCallRecord.ComplexEventInterval interval : intervals) {
            if (interval.getIntervalSVType() == type) {
                return interval;
            }
        }
        return null;
    }

    private static boolean hasSinkAttributes(final SVCallRecord record) {
        final Map<String, Object> attrs = record.getAttributes();
        return attrs.containsKey(GATKSVVCFConstants.SINK_CHROM)
                && attrs.containsKey(GATKSVVCFConstants.SINK_POS)
                && attrs.containsKey(GATKSVVCFConstants.SINK_END);
    }

    private static String getSinkChrom(final SVCallRecord record) {
        return String.valueOf(record.getAttributes().get(GATKSVVCFConstants.SINK_CHROM));
    }

    private static int getSinkPos(final SVCallRecord record) {
        final Object val = record.getAttributes().get(GATKSVVCFConstants.SINK_POS);
        return val instanceof Number ? ((Number) val).intValue() : Integer.parseInt(String.valueOf(val));
    }

    private static int getSinkEnd(final SVCallRecord record) {
        final Object val = record.getAttributes().get(GATKSVVCFConstants.SINK_END);
        return val instanceof Number ? ((Number) val).intValue() : Integer.parseInt(String.valueOf(val));
    }
}
