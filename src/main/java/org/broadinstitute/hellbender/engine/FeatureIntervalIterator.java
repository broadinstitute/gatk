package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Iterator implementation of Feature traversal by intervals.
 *
 * Given a List of intervals, queries our backing data source for Features overlapping
 * each successive interval and iterates over them, while also guaranteeing that each
 * Feature overlapping our intervals will only be returned once during the entire iteration.
 *
 * Requires that the provided List of intervals consist of non-overlapping intervals
 * sorted in increasing order of start position.
 */
class FeatureIntervalIterator<T extends Feature> implements CloseableTribbleIterator<T> {
    private final String sourceName;
    private final FeatureReader<T> featureReader;
    private final Iterator<SimpleInterval> intervalIterator;
    private CloseableTribbleIterator<T> featuresInCurrentInterval;
    private T nextFeature;
    private SimpleInterval currentInterval;
    private SimpleInterval previousInterval;

    /**
     * Initialize a FeatureIntervalIterator with a set of intervals.
     *
     * Requires that the provided List of intervals consist of non-overlapping intervals
     * sorted in increasing order of start position.
     *
     * @param intervals intervals to use for traversal. must be non-overlapping and sorted by increasing start position.
     * @param featureReader Feature reader used to perform queries
     * @param sourceName a logical name for the backing source of Features (used for error messages)
     */
    public FeatureIntervalIterator( final List<SimpleInterval> intervals, final FeatureReader<T> featureReader, final String sourceName ) {
        this.intervalIterator = intervals.iterator();
        this.featureReader = featureReader;
        this.sourceName = sourceName;
        nextFeature = loadNextNovelFeature();
    }

    @Override
    public boolean hasNext() {
        return nextFeature != null;
    }

    @Override
    public T next() {
        if ( nextFeature == null ) {
            throw new NoSuchElementException("No more Features for current interval set");
        }

        final T toReturn = nextFeature;
        nextFeature = loadNextNovelFeature();
        return toReturn;
    }

    /**
     * @return the next Feature from our data source that we HAVEN'T previously encountered,
     *         or null if we're out of Features
     */
    private T loadNextNovelFeature() {
        T candidateFeature;

        do {
            candidateFeature = loadNextFeature();

            if ( candidateFeature != null && featureIsNovel(candidateFeature) ) {
                return candidateFeature;
            }
        } while ( candidateFeature != null );

        return null;
    }

    /**
     * @return the next Feature from our data source (regardless of whether we've encountered it before or not),
     *         or null if we're out of Features
     */
    private T loadNextFeature() {
        // If we're out of Features for the current interval, repeatedly query the next interval
        // until we find one with overlapping Features. Return null if we run out of intervals.
        while ( featuresInCurrentInterval == null || ! featuresInCurrentInterval.hasNext() ) {
            if ( ! queryNextInterval() ) {
                return null;
            }
        }

        // If we reach here, we're guaranteed to have at least one Feature left to consume in the current interval.
        return featuresInCurrentInterval.next();
    }

    /**
     * Determines whether a Feature is novel (hasn't been encountered before on this iteration). A Feature
     * hasn't been encountered before if we're either on the very first query interval, or the Feature doesn't
     * overlap our previous query interval.
     *
     * @param feature Feature to test
     * @return true if we haven't seen the Feature before on this iteration, otherwise false
     */
    private boolean featureIsNovel( final T feature ) {
        return previousInterval == null || ! previousInterval.overlaps(new SimpleInterval(feature));
    }

    /**
     * Performs a query on the next interval in our interval List, and initializes all members appropriately
     * to prepare for processing the Features overlapping that interval.
     *
     * @return true if we successfully queried the next interval, false if there are no more intervals to query
     */
    private boolean queryNextInterval() {
        // Make sure to close out the query iterator for the previous interval, since Tribble only allows us
        // to have one iterator open over our FeatureReader at a time.
        if ( featuresInCurrentInterval != null ) {
            featuresInCurrentInterval.close();
            featuresInCurrentInterval = null;
        }

        if ( ! intervalIterator.hasNext() ) {
            currentInterval = previousInterval = null;
            return false;
        }

        previousInterval = currentInterval;
        currentInterval = intervalIterator.next();
        try {
            featuresInCurrentInterval = featureReader.query(currentInterval.getContig(), currentInterval.getStart(), currentInterval.getEnd());
            return true;
        }
        catch ( IOException e ) {
            throw new GATKException("Error querying " + sourceName + " over interval " + currentInterval, e);
        }
    }

    @Override
    public Iterator<T> iterator() {
        return this;
    }

    @Override
    public void close() {
        if ( featuresInCurrentInterval != null ) {
            featuresInCurrentInterval.close();
        }
    }
}
