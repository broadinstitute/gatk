package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

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
 */
public class FeatureContext {

    /**
     * FeatureManager containing backing data sources for all discovered Feature arguments
     */
    private final FeatureManager featureManager;

    /**
     * We will return Features overlapping this interval
     */
    private final GenomeLoc interval;

    /**
     * Creates a new FeatureContext given a FeatureManager and query interval
     *
     * @param featureManager FeatureManager containing backing data sources for all discovered Feature arguments
     * @param interval interval to constrain queries on this FeatureContext
     */
    public FeatureContext( final FeatureManager featureManager, final GenomeLoc interval ) {
        this.featureManager = featureManager;
        this.interval = interval;
    }

    /**
     * Gets our query interval (the interval that all Features returned by this FeatureContext overlap)
     *
     * @return query interval for this FeatureContext
     */
    public GenomeLoc getInterval() {
        return interval;
    }

    /**
     * Gets all Features from the source represented by the provided FeatureInput argument that overlap
     * this FeatureContext's query interval.
     *
     * Returned Features are not guaranteed to be in any particular order.
     *
     * @param featureDescriptor FeatureInput argument for which to fetch Features
     * @param <T> type of Feature in the data source backing the provided FeatureInput
     * @return all Features in the data source backing the provided FeatureInput that overlap
     *         this FeatureContext's query interval
     */
    public <T extends Feature> List<T> getValues( final FeatureInput<T> featureDescriptor ) {
        return featureManager.getFeatures(featureDescriptor, interval);
    }

    /**
     * Gets all Features from the source represented by the provided FeatureInput argument that overlap
     * this FeatureContext's query interval AND that start at the specified start position.
     *
     * Returned Features are not guaranteed to be in any particular order.
     *
     * @param featureDescriptor FeatureInput argument for which to fetch Features
     * @param featureStart All returned Features must start at this position, in addition to overlapping this
     *                     FeatureContext's query interval
     * @param <T> type of Feature in the data source backing the provided FeatureInput
     * @return all Features in the data source backing the provided FeatureInput that overlap
     *         this FeatureContext's query interval AND that start at the specified start position
     */
    public <T extends Feature> List<T> getValues( final FeatureInput<T> featureDescriptor, final int featureStart ) {
        return subsetToStartPosition(getValues(featureDescriptor), featureStart);
    }

    /**
     * Gets all Features from the sources represented by the provided FeatureInput arguments that overlap
     * this FeatureContext's query interval.
     *
     * Returned Features are not guaranteed to be in any particular order, or to be globally unique
     * across all sources of Features.
     *
     * @param featureDescriptors FeatureInput arguments for which to fetch Features
     * @param <T> type of Feature in the data sources backing the provided FeatureInputs
     * @return all Features in the data sources backing the provided FeatureInputs that overlap
     *         this FeatureContext's query interval
     */
    public <T extends Feature> List<T> getValues( final Collection<FeatureInput<T>> featureDescriptors ) {
        List<T> features = new ArrayList<>();
        for ( FeatureInput<T> featureSource : featureDescriptors ) {
            features.addAll(getValues(featureSource));
        }
        return features;
    }

    /**
     * Gets all Features from the sources represented by the provided FeatureInput arguments that overlap
     * this FeatureContext's query interval, AND that start at the specified start position.
     *
     * Returned Features are not guaranteed to be in any particular order, or to be globally unique
     * across all sources of Features.
     *
     * @param featureDescriptors FeatureInput arguments for which to fetch Features
     * @param featureStart All returned Features must start at this position, in addition to overlapping this
     *                     FeatureContext's query interval
     * @param <T> type of Feature in the data sources backing the provided FeatureInputs
     * @return all Features in the data sources backing the provided FeatureInputs that overlap
     *         this FeatureContext's query interval, AND that start at the specified start position
     */
    public <T extends Feature> List<T> getValues( final Collection<FeatureInput<T>> featureDescriptors, final int featureStart ) {
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
    private <T extends Feature> List<T> subsetToStartPosition( final Collection<T> features, final int start ) {
        List<T> subsettedFeatures = new ArrayList<>(features.size());
        for ( T feature : features ) {
            if ( feature.getStart() == start ) {
                subsettedFeatures.add(feature);
            }
        }
        return subsettedFeatures;
    }
}
