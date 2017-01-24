package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;

/**
 * Stores a single integer read count value for a particular {@link SimpleInterval}
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class SimpleReadCountData extends ReadCountData {

    private SimpleInterval interval;
    private int count;
    protected final static String SIMPLE_COUNT_COLUMN_NAME = "SIMPLE_COUNT";

    /**
     * Construct an empty instance of a raw read count
     *
     * @param interval corresponding target
     */
    public SimpleReadCountData(final SimpleInterval interval) {
        this.interval = Utils.nonNull(interval, "Target cannot be null. ");
        count = 0;
    }

    /**
     * Construct an instance of raw read count given target and the raw value
     *
     * @param interval corresponding target
     * @param countValue raw read count value
     */
    protected SimpleReadCountData(final SimpleInterval interval, final int countValue) {
        this.interval = Utils.nonNull(interval, "Target cannot be null. ");
        count = ParamUtils.isPositiveOrZero(countValue, "Raw read count value cannot be negative");
    }

    @Override
    public SimpleInterval getInterval() {
        return interval;
    }

    @Override
    public void updateReadCount(final GATKRead read) {
        count++;
    }

    @Override
    public void appendCountsTo(DataLine dataLine) {
        dataLine.append(count);
    }

    @Override
    public TableColumnCollection getReadCountDataColumns() {
        return new TableColumnCollection(Arrays.asList(SIMPLE_COUNT_COLUMN_NAME));
    }

    /**
     *
     * @return
     */
    public int getCount() {
        return count;
    }
}
