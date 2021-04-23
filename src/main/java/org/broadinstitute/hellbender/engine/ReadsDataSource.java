package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;
import java.util.List;

/**
 *
 * An interface for managing traversals over sources of reads.
 *
 * Two basic operations are available:
 *
 * -Iteration over all reads, optionally restricted to reads that overlap a set of intervals
 * -Targeted queries by one interval at a time
 */
public interface ReadsDataSource extends GATKDataSource<GATKRead>, AutoCloseable {

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads that overlap the given intervals,
     * and to unmapped reads if specified.
     *
     * Calls to {@link #query} are not affected by this method.
     *
     * @param intervals Our next full traversal will return reads overlapping these intervals
     * @param traverseUnmapped Our next full traversal will return unmapped reads (this affects only unmapped reads that
     *                         have no position -- unmapped reads that have the position of their mapped mates will be
     *                         included if the interval overlapping that position is included).
     */
    void setTraversalBounds(List<SimpleInterval> intervals, boolean traverseUnmapped);

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads which overlap the given intervals.
     * Calls to {@link #query} are not affected by setting these intervals.
     *
     * @param intervals Our next full traversal will return only reads overlapping these intervals
     */
    default void setTraversalBounds(List<SimpleInterval> intervals) {
        setTraversalBounds(intervals, false);
    }

    /**
     * Restricts a traversal of this data source via {@link #iterator} to only return reads that overlap the given intervals,
     * and to unmapped reads if specified.
     *
     * Calls to {@link #query} are not affected by this method.
     *
     * @param traversalParameters set of traversal parameters to control which reads get returned by the next call
     *                            to {@link #iterator}
     */
    default void setTraversalBounds(TraversalParameters traversalParameters){
        setTraversalBounds(traversalParameters.getIntervalsForTraversal(), traversalParameters.traverseUnmappedReads());
    }

    /**
     * @return true if traversals initiated via {@link #iterator} will be restricted to reads that overlap intervals
     *         as configured via {@link #setTraversalBounds}, otherwise false
     */
    boolean traversalIsBounded();

    /**
     * @return true if this datasource supports the query() operation otherwise false.
     */
    boolean isQueryableByInterval();

    /**
     * @return An iterator over just the unmapped reads with no assigned position. This operation is not affected
     *         by prior calls to {@link #setTraversalBounds}. The underlying file must be indexed.
     */
    Iterator<GATKRead> queryUnmapped();

    /**
     * Returns the SAM header for this data source.
     *
     * @return SAM header for this data source
     */
    SAMFileHeader getHeader();

    /**
     * Get the sequence dictionary for this ReadsDataSource
     *
     * @return SAMSequenceDictionary for the reads backing this datasource.
     */
    default SAMSequenceDictionary getSequenceDictionary(){
        return getHeader().getSequenceDictionary();
    }

    /**
     * @return true if this {@code ReadsDataSource} supports multiple iterations over the data
     */
    boolean supportsSerialIteration();

    /**
     * Shut down this data source permanently, closing all iterations and readers.
     */
    @Override  //Overriden here to disallow throwing checked exceptions.
    void close();
}
