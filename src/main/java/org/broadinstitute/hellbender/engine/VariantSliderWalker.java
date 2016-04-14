package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.SliderArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;

/**
 * A VariantSliderWalker is a tool that processes a single {@link LocalReadShard} window over
 * the genome/intervals at a time, with the ability to query overlapping optional sources reference data, and/or reads/features.
 *
 * The windows are constructed over each interval with a concrete window and step size, and padding overlapping
 * sources if requested.
 *
 * SlidingWindow authors must implement the {@link #apply}, {@link #defaultWindowSize}, {@link #defaultWindowStep} and {@link
 * #defaultWindowPadding} methods, and may optionally implement {@link #onTraversalStart} and/or {@link #onTraversalSuccess}.
 *
 * For control the window arguments (window-size, window-step and window-padding) exposed to the user, authors may implement
 * {@link #provideWindowSizeArgument()}, {@link #provideWindowStepArgument()} or {@link #provideWindowPaddingArgument()}.
 *
 * @author Daniel Gómez-Sánchez (magicDGS)
 */
public abstract class VariantSliderWalker extends GATKTool {
    // NOTE: using File rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public String drivingVariantFile;

    @ArgumentCollection
    protected SliderArgumentCollection sliderArgumentCollection = new SliderArgumentCollection(
            defaultWindowSize(), provideWindowSizeArgument(), defaultWindowStep(), provideWindowStepArgument(), defaultWindowPadding(), provideWindowPaddingArgument()
    );

    /**
     * Should the tool provide a window-size argument? Tools that do should override to return {@code true}.
     */
    protected boolean provideWindowSizeArgument() {
        return false;
    }

    /**
     * Should the tool provide a window-step argument? Tools that do should override to return {@code true}.
     */
    protected boolean provideWindowStepArgument() {
        return false;
    }

    /**
     * Should the tool provide a window-padding argument? Tools that do should override to return {@code true}.
     */
    protected boolean provideWindowPaddingArgument() {
        return false;
    }

    /**
     * Get the default window size for making the windows
     *
     * @return window-size (must be a positive integer)
     */
    protected abstract int defaultWindowSize();

    /**
     * Get the window step for making the windows
     *
     * @return window-step (must be a positive integer)
     */
    protected abstract int defaultWindowStep();

    /**
     * Get the length in each direction of the window to pad
     *
     * @return number of bases to pad (must be a positive integer or 0)
     */
    protected abstract int defaultWindowPadding();


    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager in GATKTool
    //we do add the driving source to the Feature manager but we do need to treat it differently and thus this field.
    private FeatureDataSource<VariantContext> drivingVariants;
    private FeatureInput<VariantContext> drivingVariantsFeatureInput;

    private List<LocalVariantShard> windows;

    @Override
    public boolean requiresFeatures() {
        return true;
    }

    @Override
    void initializeFeatures() {
        // Note: we override this method because we don't want to set feature manager to null if there are no FeatureInputs.
        // This is because we have at least 1 source of features (namely the driving dataset).
        features = new FeatureManager(this, sliderArgumentCollection.getWindowSize() + sliderArgumentCollection.getWindowPadding());
        initializeDrivingVariants();
    }

    @SuppressWarnings("unchecked")
    private void initializeDrivingVariants() {
        drivingVariantsFeatureInput = new FeatureInput<>(drivingVariantFile, "drivingVariantFile");

        //This is the data source for the driving source of variants, which uses a cache lookahead of FEATURE_CACHE_LOOKAHEAD
        drivingVariants = new FeatureDataSource<>(drivingVariantsFeatureInput, sliderArgumentCollection.getWindowSize() + sliderArgumentCollection.getWindowPadding(), VariantContext.class);

        //Add the driving datasource to the feature manager too so that it can be queried. Setting lookahead to 0 to avoid caching.
        //Note: we are disabling lookahead here because of windowed queries that need to "look behind" as well.
        features.addToFeatureSources(0, drivingVariantsFeatureInput, VariantContext.class);

        // set the intervals for traversal the driving variants
        if (hasIntervals()) {
            drivingVariants.setIntervalsForTraversal(intervalsForTraversal);
        }
    }


    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
        // first try to get the best available sequence dictionary
        final SAMSequenceDictionary dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Tool " + this.getClass().getSimpleName() + " requires some source for sequence dictionary, but none were provided");
        }
        final List<SimpleInterval> intervals = hasIntervals() ? intervalsForTraversal : IntervalUtils.getAllIntervalsForReference(dictionary);
        windows = makeVariantShards(intervals,
                sliderArgumentCollection.getWindowSize(),
                sliderArgumentCollection.getWindowStep(),
                sliderArgumentCollection.getWindowPadding(),
                dictionary);
    }

    /**
     * Shard our intervals for traversal into LocalVariantShard
     *
     * @param intervals     unmodified intervals for traversal
     * @param windowSize    the size for each window
     * @param windowStep    the step for each window
     * @param windowPadding the padding around the window
     * @return List of {@link LocalVariantShard} objects, sharded and padded as necessary
     */
    private List<LocalVariantShard> makeVariantShards(final List<SimpleInterval> intervals, final int windowSize, final int windowStep, final int windowPadding, final SAMSequenceDictionary dictionary) {
        final List<LocalVariantShard> shards = new ArrayList<>();
        for (final SimpleInterval interval : intervals) {
            shards.addAll(LocalVariantShard.divideIntervalIntoShards(interval, windowSize, windowStep, windowPadding, drivingVariants, dictionary));
        }
        return shards;
    }

    /**
     * Returns the feature input for the driving variants file.
     */
    protected final FeatureInput<VariantContext> getDrivingVariantsFeatureInput() {
        return drivingVariantsFeatureInput;
    }

    /**
     * Implementation of window-based traversal. Subclasses can override to provide their own behavior but default
     * implementation should be suitable for most uses.
     * <p>
     * The default implementation iterates over every LocalVariantShard stored in {@link #windows} and pass all the available
     * information to {@link #apply}
     */
    @Override
    public void traverse() {
        final VariantFilter filter = makeVariantFilter();
        // Since we're processing regions rather than individual variants, tell the progress meter to check the time more frequently (every 10 regions instead of every 1000 regions).
        progressMeter.setRecordsBetweenTimeChecks(10L);
        for (final LocalVariantShard window : windows) {
            // Since variants in each window are lazily fetched, we need to pass the filter to the window instead of filtering the variants directly here
            window.setVariantFilter(filter);
            apply(window,
                    new ReadsContext(reads, window.getPaddedInterval()),
                    new ReferenceContext(reference, window.getPaddedInterval()), // use the fully-padded window to fetch overlapping data
                    new FeatureContext(features, window.getPaddedInterval()));
            progressMeter.update(window.getInterval());
        }
    }

    /**
     * Process an individual window. Must be implemented by tool authors. In general, tool authors should simply stream
     * their output from apply(), and maintain as little internal state as possible.
     *
     * @param variantShard     Variants overlapping the current window (including padded regions). Will be an empty, but
     *                         non-null, context object if there is no backing source of variant data (in which case all
     *                         queries on it will return an empty array/iterator).
     * @param readsContext     Reads spanning the current window (included padded regions). Will be an empty, but non-null,
     *                         context object if there is no backing source of reads data (in which case all queries on
     *                         it will return an empty array/iterator)
     * @param referenceContext Reference bases spanning the current window (including padded regions). Will be an empty,
     *                         but non-null, context object if there is no backing source of reference data (in which
     *                         case all queries on it will return an empty array/iterator). Can request extra bases of
     *                         context around the current variant's interval by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext   Features spanning the current window (including padded regions). Will be an empty, but
     *                         non-null, context object if there is no backing source of Feature data (in which case all
     */
    protected abstract void apply(Shard<VariantContext> variantShard, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext);

    /**
     * Returns the variant filter (simple or composite) that will be applied to the variants before calling {@link #apply}.
     * The default implementation filters nothing.
     * Default implementation of {@link #traverse()} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can extend to provide own filters (ie override and call super).
     * Multiple filters can be composed by using {@link VariantFilter} composition methods.
     */
    protected VariantFilter makeVariantFilter() {
        return VariantFilterLibrary.ALLOW_ALL_VARIANTS;
    }

    /**
     * Close all data sources.
     * 
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalSuccess} and/or
     * {@link #closeTool} instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if (drivingVariants != null)
            drivingVariants.close();
    }

}
