package org.broadinstitute.hellbender.tools.gvs.common;

import htsjdk.samtools.util.Interval;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FastOverlapDetector {
    private Map<String, ArrayList<WeightedInterval>> weightsByContig = new HashMap<>();
    private int blockSize;

    public static FastOverlapDetector create(final List<WeightedInterval> intervals) {
        final FastOverlapDetector fastOverlapDetector = new FastOverlapDetector();
        // grab the first value for init purposes
        WeightedInterval startingInterval = intervals.get(0);
        // if the list is empty, something has likely gone very wrong here but still exit with an object.
        if (startingInterval == null) {
            return fastOverlapDetector;
        }

        fastOverlapDetector.blockSize = startingInterval.length();

        String currentContig = startingInterval.getContig();
        ArrayList<WeightedInterval> weightsForCurrentContig = new ArrayList<>();

        for (WeightedInterval interval : intervals) {
            if (!interval.getContig().equals(currentContig)) {
                fastOverlapDetector.weightsByContig.put(currentContig, weightsForCurrentContig);

                currentContig = interval.getContig();
                weightsForCurrentContig = new ArrayList<>();
            }
            weightsForCurrentContig.add(interval);
        }

        // clean up the last one
        fastOverlapDetector.weightsByContig.put(currentContig, weightsForCurrentContig);

        return fastOverlapDetector;
    }

    public List<WeightedInterval> getOverlaps(final Interval interval) {
        int startingBlock = interval.getStart() / blockSize;

        // There are several quirks in how intervals work.  The interval 1000-2000, when passed as an input, gets
        // interpreted as [start, end] so it would have 1001 elements.  This different from how the bed file of weights
        // is read in and interpreted, which means that an input interval of [1000, 2000] will intersect a weight
        // interval of (1000, 2000] and get the wrong answer for the first location.  The first overlapping location
        // returned will be 1001.  I work around it by detecting this boundary condition and, when possible, returning
        // the interval before the overlap that the straightforward math would otherwise give us.  This isn't a problem
        // at the very start of a chromosome because, by convention, the first location on a chromosome is position 1.
        if ((interval.getStart() % blockSize) == 0) {
            startingBlock = Math.max(0, startingBlock - 1);
        }


        // We do a -1 on the interval end because intervals are (start,end].  So one that ends on an interval boundary
        // like 140000 will not actually intersect the interal that falls in the 14000 bucket, as it goes from
        // [14001, (14000 + blockSize)].  So anyhow, subtract one from the end of the interval before placing it in a
        // bucket.  Then add 1 when you find the bucket, because otehrwise the subList command will not return the
        // correct items
        int endingBlock = (interval.getEnd() - 1) / blockSize + 1;

        return weightsByContig.get(interval.getContig()).subList(startingBlock, endingBlock);
    }
}
