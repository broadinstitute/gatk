package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Wrapper around FeatureManager that presents Feature data from a particular interval to a client tool
 * without improperly exposing engine internals.
 *
 * The client passes in one or more FeatureInputs that were declared as tool arguments, and gets back a List
 * of all Features from those FeatureInputs overlapping the interval spanned by this FeatureContext.
 *
 * Features returned may optionally be additionally constrained to start at a particular position.
 *
 * Features are returned strongly-typed based on the type parameter of each FeatureInput requested,
 * so a query on a FeatureInput<VariantContext> will return VariantContext objects (no casting
 * required by tool authors).
 *
 * Feature sources are lazily queried, so there's no overhead if the client chooses not to examine
 * the FeatureContext it's passed.
 *
 * A FeatureContext may have no backing data source and/or interval. In these cases, queries on it will always
 * return empty Lists. You can determine whether there is a backing source of Features via
 * {@link #hasBackingDataSource()}, and whether there is an interval via {@link #getInterval}
 *
 * Note: This class is NOT intended to be extended outside of the testing harness.
 */
@DoNotSubclass
public class FeatureContext {

    /**
     * FeatureManager containing backing data sources for all discovered Feature arguments.
     * Null if there are no sources of Features.
     */
    private final FeatureManager featureManager;

    /**
     * We will return Features overlapping this interval. Null if this context has no known location
     * (eg., we are dealing with unmapped data).
     */
    private final SimpleInterval interval;

    /**
     * Creates an empty FeatureContext with no backing data source. All queries on this context will
     * return an empty List.
     */
    public FeatureContext() {
        this(null, null);
    }

    /**
     * Creates a new FeatureContext given a FeatureManager and a query interval. These may be null if
     * no sources of Features are available and/or we don't have a known location on the reference,
     * in which case all queries on this context will return an empty List.
     *
     * @param featureManager FeatureManager containing backing data sources for all discovered Feature arguments. Null if there are no sources of Features.
     * @param interval Interval to constrain queries on this FeatureContext. Null if we have no known location.
     */
    public FeatureContext(final FeatureManager featureManager, final SimpleInterval interval) {
        this.featureManager = featureManager;
        this.interval = interval;
    }

    /**
     * Determines whether this FeatureContext has a backing source of Features. A FeatureContext with
     * no backing data source will always return an empty List in response to a query.
     *
     * @return true if this FeatureContext has a backing source of Features, otherwise false
     */
    public boolean hasBackingDataSource() {
        return featureManager != null;
    }

    /**
     * Gets our query interval (the interval that all Features returned by this FeatureContext overlap).
     * Null if this context has no interval.
     *
     * @return query interval for this FeatureContext (may be null)
     */
    public SimpleInterval getInterval() {
        return interval;
    }

    /**
     * Gets the header associated with the provided FeatureInput
     *
     * @param featureDescriptor FeatureInput whose header we want to retrieve
     * @param <T> type of Feature in our FeatureInput
     * @return header for the provided FeatureInput (null if we have no backing data sources)
     */
    public <T extends Feature> Object getHeader(final FeatureInput<T> featureDescriptor) {
        return featureManager != null ? featureManager.getHeader(featureDescriptor) : null;
    }

    /**
     * Gets all Features from the source represented by the provided FeatureInput argument that overlap
     * this FeatureContext's query interval. Will return an empty List if this FeatureContext has
     * no backing source of Features and/or interval.
     *
     * Returned Features are not guaranteed to be in any particular order.
     *
     * @param featureDescriptor FeatureInput argument for which to fetch Features
     * @param <T> type of Feature in the data source backing the provided FeatureInput
     * @return All Features in the data source backing the provided FeatureInput that overlap
     *         this FeatureContext's query interval. Empty List if there is no backing data source and/or interval.
     */
    public <T extends Feature> List<T> getValues(final FeatureInput<T> featureDescriptor) {
        return getValues(featureDescriptor, interval);
    }

    /**
     * Gets all Features from the source represented by the provided FeatureInput argument that overlap this
     * FeatureContext's query interval as expanded by the specified number of leading/trailing bases.
     * Returns an empty List if this FeatureContext has no backing source of Features and/or interval.
     *
     * Returned Features are not guaranteed to be in any particular order.
     *
     * Note: if windowLeadingBases > 0 and query lookahead caching is enabled for the underlying FeatureDataSource,
     * there will be a cache miss on every call, since the current caching scheme prefetches Features after
     * the current query interval but not before it (on the assumption that
     * the common access pattern involves gradually increasing query intervals).
     *
     * @param featureDescriptor FeatureInput argument for which to fetch Features
     * @param <T> type of Feature in the data source backing the provided FeatureInput
     * @param windowLeadingBases Number of extra reference bases to include before the start of our interval. Must be >= 0.
     * @param windowTrailingBases Number of extra reference bases to include after the end of our interval. Must be >= 0.
     * @return All Features in the data source backing the provided FeatureInput that overlap
     *         this FeatureContext's query interval as expanded by the specified number of leading/trailing bases.
     *         Empty List if there is no backing data source and/or interval.
     */
    public <T extends Feature> List<T> getValues(final FeatureInput<T> featureDescriptor, final int windowLeadingBases, final int windowTrailingBases) {
        return getValues(featureDescriptor, getQueryInterval(windowLeadingBases, windowTrailingBases));
    }

    /**
     * Gets all Features from the source represented by the provided FeatureInput argument that overlap the given interval.
     * Returns an empty List if this FeatureContext has no backing source of Features and/or interval.
     *
     * Returned Features are not guaranteed to be in any particular order.
     *
     * Note: if query lookahead caching is enabled for the underlying FeatureDataSource,
     * there will be a cache miss on almost every call, since the current caching scheme prefetches Features after
     * the current query interval but not before it (on the assumption that
     * the common access pattern involves gradually increasing query intervals).
     *
     * @param featureDescriptor FeatureInput argument for which to fetch Features
     * @param <T> type of Feature in the data source backing the provided FeatureInput
     * @return All Features in the data source backing the provided FeatureInput that overlap
     *         this FeatureContext's query interval as expanded by the specified number of leading/trailing bases.
     *         Empty List if there is no backing data source and/or interval.
     */
    public <T extends Feature> List<T> getValues(final FeatureInput<T> featureDescriptor, final SimpleInterval queryInterval) {
        if (featureManager == null || queryInterval == null || featureDescriptor == null) {
            return Collections.emptyList();
        }
        return featureManager.getFeatures(featureDescriptor, queryInterval);
    }

    /**
     * Gets the query interval expanded by the specified number of leading/trailing bases, or null if this context has no interval.
     *
     * @param windowLeadingBases Number of extra bases to include before the start of our interval. Must be >= 0.
     * @param windowTrailingBases Number of extra bases to include after the end of our interval. Must be >= 0.
     * @return full expanded window of bases spanned by this context as a SimpleInterval
     *         (will be null if this context has no interval)
     */
    private SimpleInterval getQueryInterval(final int windowLeadingBases, final int windowTrailingBases){
        Utils.validateArg(windowLeadingBases >= 0, "Window starts after the current interval");
        Utils.validateArg(windowTrailingBases >= 0, "Window ends before the current interval");

        if (interval == null) {
            return null;
        } else if (windowLeadingBases == 0 && windowTrailingBases == 0){
            return interval;
        }
        return new SimpleInterval(interval.getContig(), windowStart(interval, windowLeadingBases), windowStop(interval, windowTrailingBases));
    }

    /**
     * Determines the start of the expanded query window, bounded by 1.
     *
     * @param locus The locus to expand.
     * @param windowLeadingBases number of bases to attempt to expand relative to the locus start (>= 0)
     * @return The start of the expanded window.
     */
    private int windowStart(final SimpleInterval locus, final int windowLeadingBases) {
        return Math.max(locus.getStart() - windowLeadingBases, 1);
    }

    /**
     * Determines the stop of the expanded query window.
     *
     * @param locus The locus to expand.
     * @param windowTrailingBases number of bases to attempt to expand relative to the locus end (>= 0)
     * @return The end of the expanded window.
     */
    private int windowStop(final SimpleInterval locus, final int windowTrailingBases) {
        //Note: queries past the end of contig are handled but overflow
        // blows up too late so we change it here to blow up early.
        return Math.addExact(locus.getEnd(), windowTrailingBases);//blow up on overflow
    }

    /**
     * Gets all Features from the source represented by the provided FeatureInput argument that overlap
     * this FeatureContext's query interval AND that start at the specified start position.
     * Will return an empty List if this FeatureContext has no backing source of Features and/or interval.
     *
     * Returned Features are not guaranteed to be in any particular order.
     *
     * @param featureDescriptor FeatureInput argument for which to fetch Features
     * @param featureStart All returned Features must start at this position, in addition to overlapping this
     *                     FeatureContext's query interval
     * @param <T> type of Feature in the data source backing the provided FeatureInput
     * @return All Features in the data source backing the provided FeatureInput that overlap
     *         this FeatureContext's query interval AND that start at the specified start position.
     *         Empty List if there is no backing data source and/or interval.
     */
    public <T extends Feature> List<T> getValues(final FeatureInput<T> featureDescriptor, final int featureStart) {
        if (featureManager == null || interval == null) {
            return Collections.emptyList();
        }

        return subsetToStartPosition(getValues(featureDescriptor), featureStart);
    }

    /**
     * Gets all Features from the sources represented by the provided FeatureInput arguments that overlap
     * this FeatureContext's query interval. Will return an empty List if this FeatureContext has no
     * backing source of Features and/or interval.
     *
     * Returned Features are not guaranteed to be in any particular order, or to be globally unique
     * across all sources of Features.
     *
     * @param featureDescriptors FeatureInput arguments for which to fetch Features
     * @param <T> type of Feature in the data sources backing the provided FeatureInputs
     * @return All Features in the data sources backing the provided FeatureInputs that overlap
     *         this FeatureContext's query interval. Empty List if there is no backing data source and/or interval.
     */
    public <T extends Feature> List<T> getValues(final Collection<FeatureInput<T>> featureDescriptors) {
        if (featureManager == null || interval == null || featureDescriptors.isEmpty()) {
            return Collections.emptyList();
        }

        final List<T> features = new ArrayList<>();
        for (FeatureInput<T> featureSource : featureDescriptors) {
            features.addAll(getValues(featureSource));
        }
        return features;
    }

    /**
     * Gets all Features from the sources represented by the provided FeatureInput arguments that overlap
     * this FeatureContext's query interval, AND that start at the specified start position. Will return
     * an empty List if this FeatureContext has no backing source of Features and/or interval.
     *
     * Returned Features are not guaranteed to be in any particular order, or to be globally unique
     * across all sources of Features.
     *
     * @param featureDescriptors FeatureInput arguments for which to fetch Features
     * @param featureStart All returned Features must start at this position, in addition to overlapping this
     *                     FeatureContext's query interval
     * @param <T> type of Feature in the data sources backing the provided FeatureInputs
     * @return All Features in the data sources backing the provided FeatureInputs that overlap
     *         this FeatureContext's query interval, AND that start at the specified start position.
     *         Empty List if there is no backing data source and/or interval.
     */
    public <T extends Feature> List<T> getValues(final Collection<FeatureInput<T>> featureDescriptors, final int featureStart) {
        if (featureManager == null || interval == null) {
            return Collections.emptyList();
        }

        return subsetToStartPosition(getValues(featureDescriptors), featureStart);
    }

    /**
     * Helper method to subset a list of Features to only those that start at a particular start position
     *
     * @param features list of Features to subset
     * @param start the required start position for returned Features
     * @param <T> type of Feature we're dealing with
     * @return List of all Features from the features list that start at the specified start position
     */
    private <T extends Feature> List<T> subsetToStartPosition(final Collection<T> features, final int start) {
        return features.stream().filter(feat -> feat.getStart() == start).collect(Collectors.toList());
    }
}
