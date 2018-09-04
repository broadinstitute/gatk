package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;

/**
 * A FeatureWalker is a tool that processes a {@link Feature} at a time from a source of Features, with
 * optional contextual information from a reference, sets of reads, and/or supplementary sources
 * of Features.

 * Subclasses must implement the {@link #apply(Feature, ReadsContext, ReferenceContext, FeatureContext)} method to process each Feature,
 * as well as {@link #isAcceptableFeatureType(Class)} and {@link #getDrivingFeatureFile()}, and may optionally implement
 * {@link #onTraversalStart()}, {@link #onTraversalSuccess()}, and/or {@link #closeTool()}.
 *
 * @param <F> the driving feature type.
 */
public abstract class FeatureWalker<F extends Feature> extends GATKTool {

    private FeatureDataSource<F> drivingFeatures;

    @Override
    public boolean requiresFeatures(){
        return true;
    }
    
    @Override
    public String getProgressMeterRecordLabel() { return "features"; }

    @Override
    void initializeFeatures() {
        features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                      referenceArguments.getReferencePath());
        initializeDrivingFeatures();
    }

    /**
     * Set the intervals for traversal in the driving features.
     *
     * Marked final so that subclasses don't override it. Subclasses should override {@link #onTraversalStart} instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        // set the intervals for the feature here, because they are not initialized when initialize features is set
        if ( hasUserSuppliedIntervals() ) {
            drivingFeatures.setIntervalsForTraversal(userIntervals);
        }
    }

    @SuppressWarnings("unchecked")
    private void initializeDrivingFeatures() {
        final File drivingFile = getDrivingFeatureFile();
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(drivingFile);
        if (isAcceptableFeatureType(codec.getFeatureType())) {
            drivingFeatures = new FeatureDataSource<>(new FeatureInput<>(drivingFile.getAbsolutePath()), FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, null, cloudPrefetchBuffer, cloudIndexPrefetchBuffer, referenceArguments.getReferencePath());

            final FeatureInput<F> drivingFeaturesInput = new FeatureInput<>(drivingFile.getAbsolutePath(), "drivingFeatureFile");
            features.addToFeatureSources(0, drivingFeaturesInput, codec.getFeatureType(), cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                         referenceArguments.getReferencePath());
        } else {
            throw new UserException("File " + drivingFile + " contains features of the wrong type.");
        }
    }

    /**
     * Returns whether the given class of features is acceptable for this walker.
     */
    protected abstract boolean isAcceptableFeatureType(Class<? extends Feature> featureType);

    /**
     * Implementation of Feature-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     */
    @Override
    public void traverse() {
        CountingReadFilter readFilter = makeReadFilter();
        // Process each feature in the input stream.
        Utils.stream(drivingFeatures).forEach(feature -> {
                    final SimpleInterval featureInterval = new SimpleInterval(feature);
                    apply(feature,
                            new ReadsContext(reads, featureInterval, readFilter),
                            new ReferenceContext(reference, featureInterval),
                            new FeatureContext(features, featureInterval));
                    progressMeter.update(feature);
                });
    }

    /**
     * Process an individual feature.
     * In general, subclasses should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param feature Current Feature being processed.
     * @param readsContext Reads overlapping the current feature. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param referenceContext Reference bases spanning the current feature. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current feature's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current feature. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply(final F feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext );

    /**
     * Close the reads and reference data sources.
     *
     * Marked final so that subclasses don't override it. Subclasses should override {@link #onTraversalSuccess()} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( drivingFeatures != null ) {
            drivingFeatures.close();
        }
    }

    /**
     * Returns the file that contains the driving features.
     *
     * @return never {@code null}.
     */
    public abstract File getDrivingFeatureFile();
}
