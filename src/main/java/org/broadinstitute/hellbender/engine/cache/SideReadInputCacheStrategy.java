package org.broadinstitute.hellbender.engine.cache;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.function.Function;

/**
 * A {@code LocatableCacheStrategy} used to cache {@code GATKRead} objects that are side inputs. Assumes serial
 * queries will be over increasing intervals. Performance will suffer if a subsequent query interval covers territory
 * from a position earlier than a previous query.
 *
 * NOTE: this implementation uses the same criteria as BAMFileReader to determine whether to return the unmapped
 * read in a mate pair where one read is mapped and overlaps the query interval and one read is unmapped but placed
 * (the unmapped mate is only returned if the *start position* overlaps the query; otherwise only the mapped read
 * of the pair is returned). For cache trimming purposes however, the two reads will always be trimmed at the same
 * time.
 *
 * @param <CACHED_READ> Type of Locatable being cached.
 */
public class SideReadInputCacheStrategy<CACHED_READ extends GATKRead> implements LocatableCacheStrategy<CACHED_READ> {
    private static final Logger logger = LogManager.getLogger(LocatableCache.class);
    private static final int EXPECTED_MAX_OVERLAPPING_READS_DURING_CACHE_TRIM = 128;

    private int lookAheadBases;
    final Function<SimpleInterval, Iterator<CACHED_READ>> queryResultsProvider;

    /**
     * @param lookAheadBases number of bases beyond the requested interval to cache
     * @param queryResultsProvider @{code Function} that takes a query interval and returns an iterator of
     * {@code Locatable} over that interval
     */
    public SideReadInputCacheStrategy(final int lookAheadBases, final Function<SimpleInterval, Iterator<CACHED_READ>> queryResultsProvider) {
        this.lookAheadBases = lookAheadBases;
        this.queryResultsProvider = queryResultsProvider;
    }

    @Override
    public SimpleInterval getCacheIntervalFromQueryInterval(SimpleInterval queryInterval) {
        return new SimpleInterval(queryInterval.getContig(), queryInterval.getStart(), Math.addExact(queryInterval.getEnd(), lookAheadBases));
    }

    /**
     * Return reads from the cache that overlap a query interval.
     *
     * NOTE: this implementation uses the same criteria as BAMFileReader to determine whether to return the unmapped
     * read in a mate pair where one read is mapped and overlaps the query interval and one read is unmapped but placed
     * (the unmapped mate is only returned if the *start position* overlaps the query; otherwise only the mapped read
     * of the pair is returned). For cache trimming purposes however, the two reads will always be trimmed at the same
     * time.
     *
     * @param cache the cache object
     * @param queryInterval interval being queried
     * @return reads that overlap the query interval
     */
    public List<CACHED_READ> queryCache(final Deque<CACHED_READ> cache, final SimpleInterval queryInterval) {
        List<CACHED_READ> matchingReads = new ArrayList<>(cache.size());

        // Find (but do not remove from our cache) all reads that start before or on the provided stop position
        for ( CACHED_READ candidateRead : cache ) {
            if ( candidateRead.getAssignedStart() > queryInterval.getEnd() ) {
                break; // No more possible matches among the remaining cached Reads, so stop looking
            }
            if (candidateRead.isPaired() && candidateRead.isUnmapped()) {
                // in order to keep the results identical to those that would have been returned without caching,
                // use the start position of the unmapped read
                if (candidateRead.getAssignedStart() >= queryInterval.getStart()) {
                    matchingReads.add(candidateRead);
                }
            } else {
                matchingReads.add(candidateRead);
            }
        }
        return matchingReads;
    }

    @Override
    public Iterator<CACHED_READ> refillCache(SimpleInterval queryInterval) {
        return queryResultsProvider.apply(queryInterval);
    }

    @Override
    public SimpleInterval trimCache(final Deque<CACHED_READ> cache, final SimpleInterval cachedInterval, final SimpleInterval interval) {
        if ( interval.getStart() > cachedInterval.getEnd() ) {
            throw new GATKException(String.format("BUG: attempted to trimCache cache to an improper new start position (%d). Cache stop = %d",
                    interval.getStart(), cachedInterval.getEnd()));
        }

        List<CACHED_READ> overlappingReadsBeforeNewStart = new ArrayList<>(EXPECTED_MAX_OVERLAPPING_READS_DURING_CACHE_TRIM);

        // In order to trim the cache to the new start position, we need to find all reads in the cache that start
        // before the new start position, and discard those that don't overlap the new start while keeping those
        // that do overlap. We can stop once we find a read that starts on or after the new start position, since
        // the reads are assumed to be sorted by start position.
        //
        // For mate pairs where one read is unmapped, we need to keep the pairs together, so we use the territory
        // covered by the mapped mate to determine whether or not to trim the pair.
        while ( ! cache.isEmpty() && cache.getFirst().getAssignedStart() < interval.getStart() ) {
            CACHED_READ readBeforeNewStart = cache.removeFirst();
            CACHED_READ matedRead = null;
            if (readBeforeNewStart.isUnmapped()) {
                // we found an unmapped read who's mate has not been seen yet
                matedRead = findMappedMateForUnmappedRead(readBeforeNewStart, cache);
            }
            else if (readBeforeNewStart.isPaired() && readBeforeNewStart.mateIsUnmapped()) {
                // we found a mapped read that has an unmapped mate
                matedRead = findUnmappedMateForMappedRead(readBeforeNewStart, cache);
            }
            // Our trim criteria for a pair, one of which is unmapped, should be based on the alignment end from the
            // MAPPED mate of the pair if there is one. Since the mates can appear in either order, determine which
            // read is the mapped mapped and use it's value for the end calculation.
            int mappedEnd = readBeforeNewStart.isUnmapped() && matedRead != null ?
                    matedRead.getEnd() :
                    readBeforeNewStart.getEnd();
            if ( mappedEnd >= interval.getStart() ) {
                overlappingReadsBeforeNewStart.add(readBeforeNewStart);
                if (matedRead != null)
                    // mated reads should travel together
                    overlappingReadsBeforeNewStart.add(matedRead);
            }
        }

        // Add back the reads that started before the new start but overlapped it in the reverse of the order in
        // which we encountered them so that their original relative ordering in the cache is restored.
        for ( int i = overlappingReadsBeforeNewStart.size() - 1; i >= 0; --i ) {
            cache.addFirst(overlappingReadsBeforeNewStart.get(i));
        }

        // Record our new start boundary
        return new SimpleInterval(cachedInterval.getContig(), interval.getStart(), cachedInterval.getEnd());
    }

    // Find the unmapped read's mapped mate (if its next in the cache)
    private CACHED_READ findMappedMateForUnmappedRead(
            final GATKRead unmappedRead,
            final Deque<CACHED_READ> cache) {
        if (!unmappedRead.isPaired()) {
            // We should never find an unmapped read in the cache that is not paired
            // TODO: is throwing too strict here ?
            throw new GATKException(String.format("Found unmapped, unpaired read: '%s' with no mate", unmappedRead));
        }
        if (!cache.isEmpty()) {
            final GATKRead nextRead = cache.getFirst();
            if (nextRead.isPaired() && !nextRead.isUnmapped() && nextRead.getName().equals(unmappedRead.getName())) {
                return cache.removeFirst();
            }
        }
        // Its possible that the input source actually doesn't contain the mated read, which we tolerate
        logger.info(String.format("An unmapped read '%s' with no corresponding, mapped, paired mate was found", unmappedRead));
        return null;
    }

    // Find the read's unmapped mate (if its next in the cache).
    private CACHED_READ findUnmappedMateForMappedRead(
            final GATKRead readWithUnmappedMate,
            final Deque<CACHED_READ> cache)
    {
        if (cache.isEmpty()) {
            logger.info(String.format("The unmapped mate of mapped read '%s' is missing from the input", readWithUnmappedMate));
        } else {
            final GATKRead nextRead = cache.getFirst();
            if (!nextRead.isUnmapped() || !nextRead.getName().equals(readWithUnmappedMate.getName())) {
                logger.info(String.format("A mapped read '%s' with no corresponding paired, unmapped mate was found", readWithUnmappedMate));
            }
            else {
                return cache.removeFirst();
            }
        }
        return null;
    }

}
