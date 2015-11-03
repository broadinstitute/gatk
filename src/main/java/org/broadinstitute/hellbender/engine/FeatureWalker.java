package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.Collections;
import java.util.stream.StreamSupport;

/**
 * General Feature walker.
 * <p>
 *     For tools that iterate along a list of features.
 * </p>
 *
 * @param <F> the driving feature type.
 */
public abstract class FeatureWalker<F extends Feature> extends GATKTool {

    private FeatureDataSource<F> drivingFeatures;

    @Override
    void initializeFeatures() {
        features = new FeatureManager(this);
        initializeDrivingFeatures();
    }

    @SuppressWarnings("unchecked")
    private void initializeDrivingFeatures() {
        final File drivingFile = getDrivingFeatureFile();
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(drivingFile);
        if (isAcceptableFeatureType(codec.getFeatureType())) {
            drivingFeatures = new FeatureDataSource<>(drivingFile, (FeatureCodec<F, ?>)codec);

            final FeatureInput<F> drivingFeaturesInput = new FeatureInput<>("drivingVariantFile", Collections.emptyMap(), drivingFile);
            features.addToFeatureSources(0, drivingFeaturesInput, StandardArgumentDefinitions.VARIANT_LONG_NAME, StandardArgumentDefinitions.VARIANT_SHORT_NAME, codec.getFeatureType());
        } else {
            throw new UserException("File " + drivingFile + " contains features of the wrong type.");
        }

        if ( hasIntervals() ) {
            drivingFeatures.setIntervalsForTraversal(intervalsForTraversal);
        }
    }

    protected abstract boolean isAcceptableFeatureType(Class<? extends Feature> featureType);

    /**
     * Implementation of variant-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     */
    @Override
    public void traverse() {
        // Process each variant in the input stream.
        StreamSupport.stream(drivingFeatures.spliterator(), false)
                .forEach(feature -> {
                    final SimpleInterval variantInterval = new SimpleInterval(feature);
                    apply(feature,
                            new ReadsContext(reads, variantInterval),
                            new ReferenceContext(reference, variantInterval),
                            new FeatureContext(features, variantInterval));
                });
    }

    /**
     * Process an individual variant. Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param feature Current variant being processed.
     * @param readsContext Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param referenceContext Reference bases spanning the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current read's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current variant. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply(final F feature, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext );

    /**
     * Close the reads and reference data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalDone() instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( drivingFeatures != null )
            drivingFeatures.close();
    }

    /**
     * Returns the file that contains the driving locatables.
     *
     * @return never {@code null}.
     */
    public abstract File getDrivingFeatureFile();
}
