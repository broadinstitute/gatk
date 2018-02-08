package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * A container that stores and updates read count data pertaining to a single interval
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public abstract class ReadCountData implements Locatable {

    /**
     * Updates information about the interval's read count data given a new read
     * @param read GATK read
     */
    public abstract void updateReadCount(final GATKRead read);

    /**
     * Appends this instance read counts to a data-line object.
     *
     * @param dataLine the destination data-line.
     * @throws IllegalArgumentException if {@code dataLine is {@code null}}.
     * @throws IllegalStateException    if there is not enough room in {@code dataLine} from its current appending position
     *                                  to add all the record counts.
     */
    public abstract void appendCountsTo(final DataLine dataLine);

    /**
     * Get the interval associated with this record
     *
     * @return never {@code null}.
     */
    public abstract SimpleInterval getInterval();

    /**
     * Get the columns of the read count record.
     * Classes that extend {@link ReadCountData} should implement a concrete version of this method
     *
     * @return never {@code null}.
     */
    public abstract TableColumnCollection getReadCountDataColumns();

    /**
     * Return the aggregate read count associated with this record
     *
     * @return total read count
     */
    public abstract int getTotalReadCount();

    @Override
    public String getContig() {
        return getInterval().getContig();
    }

    @Override
    public int getStart() {
        return getInterval().getStart();
    }

    @Override
    public int getEnd() {
        return getInterval().getEnd();
    }
}