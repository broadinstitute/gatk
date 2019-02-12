package org.broadinstitute.hellbender.utils.activityprofile;

import com.google.common.base.Function;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import org.broadinstitute.hellbender.engine.MultiIntervalShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nullable;
import java.util.Iterator;
import java.util.List;

/**
 * An efficient way of representing a set of {@link ActivityProfileState}s in an interval; used by Spark.
 */
public class ActivityProfileStateRange {
    private final SimpleInterval interval;
    private final double[] activeProb;
    private final ActivityProfileState.Type[] resultState;
    private final double[] resultValue; // don't store as a Number since it uses more memory

    public ActivityProfileStateRange(MultiIntervalShard<?> shard, Iterator<ActivityProfileState> activityProfileStateIterator) {
        List<SimpleInterval> intervals = shard.getIntervals();
        this.interval = Iterables.getOnlyElement(intervals);
        int size = interval.size();
        this.activeProb = new double[size];
        this.resultState = new ActivityProfileState.Type[size];
        this.resultValue = new double[size];

        int i = 0;
        ActivityProfileState prev = null;
        while (activityProfileStateIterator.hasNext()) {
            ActivityProfileState next = activityProfileStateIterator.next();
            if (prev != null) {
                Utils.validate(next.getLoc().getContig().equals(prev.getLoc().getContig()), "Contigs differ");
                Utils.validate(next.getLoc().getStart() == prev.getLoc().getStart() + 1, "Out of order");
            }
            activeProb[i] = next.isActiveProb();
            resultState[i] = next.getResultState();
            // store null result value as a negative number, since negative numbers are illegal in ActivityProfileState
            resultValue[i] = next.getResultValue() == null ? Double.NEGATIVE_INFINITY : next.getResultValue().doubleValue();
            i++;
            prev = next;
        }
        Utils.validate(i == size, "Size is wrong");
    }

    public String getContig() {
        return interval.getContig();
    }

    private Iterator<ActivityProfileState> toIteratorActivityProfileState() {
        return new AbstractIterator<ActivityProfileState>() {
            int i = 0;
            @Override
            protected ActivityProfileState computeNext() {
                if (i == interval.size()) {
                    return endOfData();
                }
                int pos = interval.getStart() + i;
                double v = resultValue[i];
                ActivityProfileState state = new ActivityProfileState(new SimpleInterval(interval.getContig(), pos, pos), activeProb[i], resultState[i], v == Double.NEGATIVE_INFINITY ? null : v);
                i++;
                return state;
            }
        };
    }

    public static Iterator<ActivityProfileState> toIteratorActivityProfileState(Iterator<ActivityProfileStateRange> it) {
        Iterator<Iterator<ActivityProfileState>> iteratorOfIterators = Iterators.transform(it, new Function<ActivityProfileStateRange, Iterator<ActivityProfileState>>() {
            @Nullable
            @Override
            public Iterator<ActivityProfileState> apply(@Nullable ActivityProfileStateRange input) {
                return input.toIteratorActivityProfileState();
            }
        });
        return Iterators.concat(iteratorOfIterators);
    }
}
