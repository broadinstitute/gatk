package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.exceptions.GATKException;
import picard.util.IntervalList.IntervalListScatterer;

import java.util.List;

public class IntervalListScattererWithSubdivisionByWeight implements IntervalListScatterer {

    private static WeightedInterval assertWeighted(Interval interval){
        if(interval instanceof WeightedInterval){
            return (WeightedInterval) interval;
        } else {
            throw new GATKException("IntervalListScattererWithSubdivisionByWeight is only a valid scattering mode when using" +
                    "WeightedIntervals, please report this to the developers.");
        }
    }

    @Override
    public long intervalWeight(final Interval interval) {
        return Math.max(1, (long)(assertWeighted(interval).getWeight()));
    }

    @Override
    public long listWeight(final IntervalList intervalList) {
        return intervalList.getIntervals().stream()
                .map(IntervalListScattererWithSubdivisionByWeight::assertWeighted)
                .mapToLong((interval -> (long)interval.getWeight()))
                .sum();
    }

    @Override
    public int deduceIdealSplitWeight(final IntervalList intervalList, final int nCount) {
        return Math.max(1, (int) Math.floorDiv(listWeight(intervalList), nCount));
    }

    @Override
    public IntervalList preprocessIntervalList(IntervalList inputList) {
        return inputList.sorted();
    }

    @Override
    public List<Interval> takeSome(final Interval interval, final long idealSplitWeight, final long currentSize, final double projectSizeOfRemaining) {
        final long targetWeight = idealSplitWeight - currentSize;
        final WeightedInterval weightedInterval = assertWeighted(interval);
        if (targetWeight >= weightedInterval.getWeight()) {
            return CollectionUtil.makeList(interval, null);
        }

        // if there's 1 base left and it's heavier than allowed it gets it's own interval
        if (targetWeight < weightedInterval.getWeight() && weightedInterval.length() == 1){
            return CollectionUtil.makeList(interval, null);
        }

        if (targetWeight == 0) {
            return CollectionUtil.makeList(null, interval);
        }

        final int targetBases = Math.max(1, (int)((double)targetWeight / weightedInterval.getPerBaseWeight()));
        final int splitPoint = interval.getStart() + targetBases;


        final Interval left = new WeightedInterval(
                interval.getContig(),
                interval.getStart(),
                splitPoint - 1,
                interval.isNegativeStrand(),
                interval.getName(),
                (int)((double)targetBases * weightedInterval.getPerBaseWeight())
        );
        final Interval right = new WeightedInterval(
                interval.getContig(),
                splitPoint,
                interval.getEnd(),
                interval.isNegativeStrand(),
                interval.getName(),
                (int)((double)(weightedInterval.length() - targetBases)* weightedInterval.getPerBaseWeight())
        );
        return CollectionUtil.makeList(left, right);
    }
}
