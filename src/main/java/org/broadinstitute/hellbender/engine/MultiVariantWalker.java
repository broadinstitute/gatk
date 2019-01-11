package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.Spliterator;

/**
 * A MultiVariantWalker is a tool that processes one variant at a time, in position order, from multiple sources of
 * variants, with optional contextual information from a reference, sets of reads, and/or supplementary sources of
 * Features.
 *
 * VariantWalker authors must implement the {@link #apply} method to process each variant, and may optionally implement
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class MultiVariantWalker extends VariantWalkerBase {

    @ArgumentCollection
    protected MultiVariantInputArgumentCollection multiVariantInputArgumentCollection = getMultiVariantInputArgumentCollection();

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager
    // in GATKTool we do add the driving source to the Feature manager but we do need to treat it differently and thus this
    // field.
    private MultiVariantDataSource drivingVariants;
    private List<FeatureInput<VariantContext>> drivingVariantsFeatureInputs = new ArrayList<>(2);

    @Override
    protected SAMSequenceDictionary getSequenceDictionaryForDrivingVariants() { return drivingVariants.getSequenceDictionary(); }

    @Override
    protected Spliterator<VariantContext> getSpliteratorForDrivingVariants() { return drivingVariants.spliterator(); }

    /**
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalStart} instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();

        if ( hasUserSuppliedIntervals() ) {
            drivingVariants.setIntervalsForTraversal(userIntervals);
        }
    }

    /**
     * Return an argument collection that provides the driving variants.  This allows subclasses to override and use a different
     * argument pattern besides the default -V
     */
    protected MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MultiVariantInputArgumentCollection.DefaultMultiVariantInputArgumentCollection();
    }

    @Override
    protected void initializeDrivingVariants() {
        multiVariantInputArgumentCollection.getDrivingVariantPaths().stream().forEach(
                f -> {
                    FeatureInput<VariantContext> featureInput = new FeatureInput<>(f);
                    if (drivingVariantsFeatureInputs.contains(featureInput)) {
                        throw new IllegalArgumentException("Feature inputs must be unique: " + featureInput.toString());
                    }
                    drivingVariantsFeatureInputs.add(featureInput);

                    //Add the driving datasource to the feature manager too so that it can be queried. Setting lookahead to 0 to avoid caching.
                    //Note: we are disabling lookahead here because of windowed queries that need to "look behind" as well.
                    features.addToFeatureSources(0, featureInput, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                                 referenceArguments.getReferencePath());
                }
        );
        drivingVariants = new MultiVariantDataSource(drivingVariantsFeatureInputs, VariantWalkerBase.FEATURE_CACHE_LOOKAHEAD, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                                     referenceArguments.getReferencePath());

        //Note: the intervals for the driving variants are set in onStartup
    }

    /**
     * Returns a list of feature inputs used for the driving variants for this source.
     */
    protected final List<FeatureInput<VariantContext>> getDrivingVariantsFeatureInputs() {
        return drivingVariantsFeatureInputs;
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    @Override
    public final VCFHeader getHeaderForVariants() {
        return drivingVariants.getHeader();
    }

    /**
     * Implementation of variant-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     */
    @Override
    public void traverse() {
        final CountingReadFilter readFilter = makeReadFilter();
        // Process each variant in the input stream.
        getTransformedVariantStream( makeVariantFilter() )
                .forEach(variant -> {
                    final SimpleInterval variantInterval = new SimpleInterval(variant);
                    apply(variant,
                            new ReadsContext(reads, variantInterval, readFilter),
                            new ReferenceContext(reference, variantInterval),
                            new FeatureContext(features, variantInterval));

                    progressMeter.update(variantInterval);
                });
    }

    /**
     * Process an individual variant. Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param variant Current variant being processed.
     * @param readsContext Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                     if there is no backing source of reads data (in which case all queries on it will return
     *                     an empty array/iterator)
     * @param referenceContext Reference bases spanning the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current variant's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current variant. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext );

    /**
     * Close all data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalSuccess} and/or
     * {@link #closeTool} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();
        if (drivingVariants != null) {
            drivingVariants.close();
        }
    }

    /**
     * Get in name-sorted order a list of samples merged from the driving variants files.
     *
     * @return SampleSet merged requiring unique names from the drivingVariants
     */
    public final SortedSet<String> getSamplesForVariants() {
        return drivingVariants.getSamples();
    }
}
