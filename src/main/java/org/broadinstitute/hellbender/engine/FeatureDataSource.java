package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.*;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.utils.IndexUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Enables traversals and queries over sources of Features, which are metadata associated with a location
 * on the genome in a format supported by our file parsing framework, Tribble. Examples of Features are
 * VCF records and hapmap records.
 *
 * Two basic operations are available on this data source:
 *
 * -Iteration over all Features in this data source, optionally restricted to Features overlapping
 *  a set of intervals if intervals are provided via {@link #setIntervalsForTraversal(List)}. Traversal
 *  by a set of intervals requires the file to have been indexed using the bundled tool IndexFeatureFile.
 *  The set of intervals provided MUST be non-overlapping and sorted in increasing order of start position.
 *
 * -Targeted queries by one interval at a time. This also requires the file to have been indexed using
 *  the bundled tool IndexFeatureFile. Targeted queries by one interval at a time are unaffected by
 *  any intervals for full traversal set via {@link #setIntervalsForTraversal(List)}.
 *
 * To improve performance in the case of targeted queries by one interval at a time, this class uses a caching
 * scheme that is optimized for the common access pattern of multiple separate queries over intervals with
 * gradually increasing start positions. It optimizes for this use case by pre-fetching records immediately
 * following each interval during a query and caching them. Performance will suffer if the access pattern is
 * random, involves queries over intervals with DECREASING start positions instead of INCREASING start positions,
 * or involves lots of very large jumps forward on the genome or lots of contig switches. Query caching
 * can be disabled, if desired.
 *
 * @param <T> The type of Feature returned by this data source
 */
public final class FeatureDataSource<T extends Feature> implements GATKDataSource<T>, AutoCloseable {
    private static final Logger logger = LogManager.getLogger(FeatureDataSource.class);

    /**
     * File backing this data source. Used mainly for error messages.
     */
    private final File featureFile;

    /**
     * Feature reader used to retrieve records from our file
     */
    private final AbstractFeatureReader<T, ?> featureReader;

    /**
     * Tribble codec used by our reader to decode the records in our file
     */
    private final FeatureCodec<T, ?> codec;

    /**
     * Iterator representing an open traversal over this data source initiated via a call to {@link #iterator}
     * (null if there is no open traversal). We need this to ensure that each iterator is properly closed,
     * and to enforce the constraint (required by Tribble) that we never have more than one iterator open
     * over our feature reader.
     */
    private CloseableTribbleIterator<T> currentIterator;

    /**
     * Our intervals for traversal. If set, restricts full traversals initiated via {@link #iterator} to
     * return only Features overlapping this set of intervals. Does not affect individual queries
     * initiated via {@link #query(SimpleInterval)} and/or {@link #queryAndPrefetch(SimpleInterval)}.
     */
    private List<SimpleInterval> intervalsForTraversal;

    /**
     * Cache containing Features from recent queries initiated via {@link #query(SimpleInterval)} and/or
     * {@link #queryAndPrefetch(SimpleInterval)}. This is guaranteed to start at the start position of the
     * most recent query, but will typically end well after the end of the most recent query. Designed to
     * improve performance of the common access pattern involving multiple queries across nearby intervals
     * with gradually increasing start positions.
     */
    private final FeatureCache<T> queryCache;

    /**
     * When we experience a cache miss (ie., a query interval not fully contained within our cache) and need
     * to re-populate the Feature cache from disk to satisfy a query, this controls the number of extra bases
     * AFTER the end of our interval to fetch. Should be sufficiently large so that typically a significant number
     * of subsequent queries will be cache hits (ie., query intervals fully contained within our cache) before
     * we have another cache miss and need to go to disk again.
     */
    private final int queryLookaheadBases;

    /**
     * An (optional) logical name assigned to this data source. May be null.
     */
    private final String name;

    /**
     * True if the file backing this data source has an accompanying index file, false if it doesn't
     */
    private final boolean hasIndex;

    /**
     * Default value for queryLookaheadBases, if none is specified. This is designed to be large enough
     * so that in typical usage (ie., query intervals with gradually increasing start locations) there will
     * be a substantial number of cache hits between cache misses, reducing the number of times we need to
     * repopulate the cache from disk.
     */
    public static final int DEFAULT_QUERY_LOOKAHEAD_BASES = 1000;

    /**
     * Creates a FeatureDataSource backed by the provided file that uses the given codec to decode records
     * from that file. The data source will have no name, and will look ahead the default number of bases
     * ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES}) during queries that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param codec codec with which to decode the records from featureFile
     */
    public FeatureDataSource( final File featureFile, final FeatureCodec<T, ?> codec ) {
        this(featureFile, codec, null, DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    /**
     * Creates a FeatureDataSource backed by the provided File that uses the provided codec to decode records
     * from that file, and assigns this data source a logical name. We will look ahead the default number of bases
     * ({@link #DEFAULT_QUERY_LOOKAHEAD_BASES}) during queries that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param codec codec with which to decode the records from featureFile
     * @param name logical name for this data source (may be null)
     */
    public FeatureDataSource( final File featureFile, final FeatureCodec<T, ?> codec, final String name ) {
        this(featureFile, codec, name, DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    /**
     * Creates a FeatureDataSource backed by the provided File that uses the provided codec to decode records
     * from that file, and assigns this data source a logical name. We will look ahead the specified number of bases
     * during queries that produce cache misses.
     *
     * @param featureFile file containing Features
     * @param codec codec with which to decode the records from featureFile
     * @param name logical name for this data source (may be null)
     * @param queryLookaheadBases look ahead this many bases during queries that produce cache misses
     */
    public FeatureDataSource( final File featureFile, final FeatureCodec<T, ?> codec, final String name, final int queryLookaheadBases ) {
        if ( featureFile == null || codec == null ) {
            throw new IllegalArgumentException("FeatureDataSource cannot be created from null file/codec");
        }
        if ( queryLookaheadBases < 0 ) {
            throw new IllegalArgumentException("Query lookahead bases must be >= 0");
        }
        if ( ! featureFile.canRead() || featureFile.isDirectory() ) {
            throw new UserException.CouldNotReadInputFile("File " + featureFile.getAbsolutePath() + " does not exist, is unreadable, or is a directory");
        }

        this.featureFile = featureFile;

        try {
            // Instruct the reader factory to not require an index. We will require one ourselves as soon as
            // a query by interval is attempted.
            this.featureReader = AbstractFeatureReader.getFeatureReader(featureFile.getAbsolutePath(), codec, false);
        }
        catch ( TribbleException e ) {
            throw new GATKException("Error initializing feature reader for file " + featureFile.getAbsolutePath(), e);
        }

        this.currentIterator = null;
        this.intervalsForTraversal = null;
        this.queryCache = new FeatureCache<>();
        this.queryLookaheadBases = queryLookaheadBases;
        this.codec = codec;
        this.name = name;
        this.hasIndex = featureReader.hasIndex(); // Cache this result, as it's fairly expensive to determine
    }

    /**
     * Returns the sequence dictionary for this source of Features.
     * Uses the dictionary from the VCF header (if present) for variant inputs,
     * otherwise attempts to create a sequence dictionary from the index file (if present).
     * Returns null if no dictionary could be created from either the header or the index.
     */
    public SAMSequenceDictionary getSequenceDictionary() {
        SAMSequenceDictionary dict = null;
        final Object header = getHeader();
        if (header instanceof VCFHeader) {
            dict = ((VCFHeader) header).getSequenceDictionary();
        }
        if (dict != null && !dict.isEmpty()) {
            return dict;
        }
        if (hasIndex) {
            return IndexUtils.createSequenceDictionaryFromFeatureIndex(featureFile);
        }
        return null;
    }

    /**
     * Restricts traversals of this data source via {@link #iterator} to only return Features that overlap the provided
     * intervals. Calls to {@link #query(SimpleInterval)} and/or {@link #queryAndPrefetch(SimpleInterval)} are not
     * affected by these intervals.
     *
     * Intervals MUST be non-overlapping and sorted in order of increasing start position, otherwise traversal
     * results will be incorrect.
     *
     * Passing in a null or empty interval List clears the intervals for traversal, making future iterations
     * over this data source unrestricted by intervals.
     *
     * @param intervals Our next full traversal will return only Features overlapping these intervals
     */
    public void setIntervalsForTraversal( final List<SimpleInterval> intervals ) {
        // Treat null and empty interval lists the same
        intervalsForTraversal = (intervals != null && !intervals.isEmpty()) ? intervals : null;

        if ( intervalsForTraversal != null && ! hasIndex ) {
            throw new UserException("File " + featureFile.getAbsolutePath() + " requires an index to enable traversal by intervals. " +
                                    "Please index this file using the bundled tool " + IndexFeatureFile.class.getSimpleName());
        }
    }


    /**
     * Gets an iterator over all Features in this data source, restricting traversal to Features
     * overlapping our intervals if intervals were provided via {@link #setIntervalsForTraversal(List)}
     *
     * Calling this method invalidates (closes) any previous iterator obtained from this method.
     *
     * @return an iterator over all Features in this data source, limited to Features that overlap the intervals supplied via {@link #setIntervalsForTraversal(List)} (if intervals were provided)
     */
    @Override
    public Iterator<T> iterator() {
        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        try {
            // Save the iterator returned so that we can close it properly later
            currentIterator = intervalsForTraversal != null ? new FeatureIntervalIterator<T>(intervalsForTraversal, featureReader, featureFile.getAbsolutePath())
                                                            : featureReader.iterator();
            return currentIterator;
        }
        catch ( IOException e ) {
            throw new GATKException("Error creating iterator over file " + featureFile.getAbsolutePath(), e);
        }
    }

    /**
     * Gets an iterator over all Features in this data source that overlap the provided interval.
     *
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     *
     * Requires the backing file to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     *
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     *
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Features overlapping this interval
     * @return an iterator over all Features in this data source that overlap the provided interval
     */
    @Override
    public Iterator<T> query( final SimpleInterval interval ) {
        return queryAndPrefetch(interval).iterator();
    }

    /**
     * Returns a List of all Features in this data source that overlap the provided interval.
     *
     * This operation is not affected by intervals provided via {@link #setIntervalsForTraversal(List)}.
     *
     * Requires the backing file to have been indexed using the IndexFeatureFile tool, and to
     * be sorted in increasing order of start position for each contig.
     *
     * Query results are cached to improve the performance of future queries during typical access
     * patterns. See notes to the class as a whole for a description of the caching strategy.
     *
     * Calling this method potentially invalidates (closes) any other open iterator obtained
     * from this data source via a call to {@link #iterator}
     *
     * @param interval retrieve all Features overlapping this interval
     * @return a List of all Features in this data source that overlap the provided interval
     */
    public List<T> queryAndPrefetch( final SimpleInterval interval ) {
        if ( ! hasIndex ) {
            throw new UserException("File " + featureFile.getAbsolutePath() + " requires an index to enable queries by interval. " +
                                    "Please index this file using the bundled tool " + IndexFeatureFile.class.getSimpleName());
        }

        // If the query can be satisfied using existing cache contents, prepare for retrieval
        // by discarding all Features at the beginning of the cache that end before the start
        // of our query interval.
        if ( queryCache.cacheHit(interval) ) {
            queryCache.trimToNewStartPosition(interval.getStart());
        }
        // Otherwise, we have a cache miss, so go to disk to refill our cache.
        else {
            refillQueryCache(interval);
        }

        // Return the subset of our cache that overlaps our query interval
        return queryCache.getCachedFeaturesUpToStopPosition(interval.getEnd());
    }

    /**
     * Refill our cache from disk after a cache miss. Will prefetch Features overlapping an additional
     * queryLookaheadBases bases after the end of the provided interval, in addition to those overlapping
     * the interval itself.
     *
     * Calling this has the side effect of invalidating (closing) any currently-open iteration over
     * this data source.
     *
     * @param interval the query interval that produced a cache miss
     */
    private void refillQueryCache( final SimpleInterval interval ) {
        // Tribble documentation states that having multiple iterators open simultaneously over the same FeatureReader
        // results in undefined behavior
        closeOpenIterationIfNecessary();

        // Expand the end of our query by the configured number of bases, in anticipation of probable future
        // queries with slightly larger start/stop positions.
        //
        // Note that it doesn't matter if we go off the end of the contig in the process, since
        // our reader's query operation is not aware of (and does not care about) contig boundaries.
        // Note: we use addExact to blow up on overflow rather than propagate negative results downstream
        final SimpleInterval queryInterval = new SimpleInterval(interval.getContig(), interval.getStart(), Math.addExact(interval.getEnd(), queryLookaheadBases));

        // Query iterator over our reader will be immediately closed after re-populating our cache
        try ( CloseableTribbleIterator<T> queryIter = featureReader.query(queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd()) ) {
            queryCache.fill(queryIter, queryInterval);
        }
        catch ( IOException e ) {
            throw new GATKException("Error querying file " + featureFile.getAbsolutePath() + " over interval " + interval, e);
        }
    }

    /**
     * Get the class of the codec being used to decode records from our file
     *
     * @return the class of the codec being used to decode records from our file
     */
    @SuppressWarnings("rawtypes")
    public final Class<? extends FeatureCodec> getCodecClass() {
        return codec.getClass();
    }

    /**
     * Get the type of Feature record stored in this data source
     *
     * @return the type of Feature record stored in this data source
     */
    public final Class<T> getFeatureType() {
        return codec.getFeatureType();
    }

    /**
     * Get the logical name of this data source. Will be null if the data source was not assigned a name.
     *
     * @return the logical name of this data source (may be null)
     */
    public String getName() {
        return name;
    }

    /**
     * Gets the header associated with this data source
     *
     * @return header associated with this data source as an Object
     */
    public Object getHeader() {
        return featureReader.getHeader();
    }

    /**
     * Permanently close this data source, invalidating any open iteration over it, and making it invalid for future
     * iterations and queries.
     */
    @Override
    public void close() {
        closeOpenIterationIfNecessary();

        logger.debug(String.format("Cache statistics for FeatureInput %s:", name != null ? name + " (" + featureFile.getName() + ")" : featureFile.getName()));
        queryCache.printCacheStatistics();

        try {
            if ( featureReader != null )
                featureReader.close();
        }
        catch ( IOException e ) {
            throw new GATKException("Error closing Feature reader for file " + featureFile.getAbsolutePath());
        }
    }

    /**
     * Close the iterator currently open over this data source, if there is one.
     */
    private void closeOpenIterationIfNecessary() {
        if ( currentIterator != null ) {
            currentIterator.close();
            currentIterator = null;
        }
    }
}
