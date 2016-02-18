package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.util.Collections;
import java.util.stream.StreamSupport;

/**
 * A VariantWalker is a tool that processes a variant at a time from a source of variants, with
 * optional contextual information from a reference, sets of reads, and/or supplementary sources
 * of Features.
 *
 * VariantWalker authors must implement the apply() method to process each read, and may optionally implement
 * onTraversalStart() and/or onTraversalSuccess().
 */
public abstract class VariantWalker extends GATKTool {

    // NOTE: using File rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public File drivingVariantFile;

    /**
     * This number controls the size of the cache for our primary and auxiliary FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 100_000;

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager in GATKTool
    //we do add the driving source to the Feature manager but we do need to treat it differently and thus this field.
    private FeatureDataSource<VariantContext> drivingVariants;
    private FeatureInput<VariantContext> drivingVariantsFeatureInput;

    /*
     * TODO: It's awkward that this traversal type requires variants yet can't override requiresFeatures() from
     * TODO: GATKTool, since the required source of variants is managed separately (as drivingVariants) from
     * TODO: our main FeatureManager in GATKTool. May need a way to register additional data sources with GATKTool.
     */

    @Override
    void initializeFeatures() {
        //Note: we override this method because we don't want to set feature manager to null if there are no FeatureInputs.
        //This is because we have at least 1 source of features (namely the driving dataset).
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD);
        initializeDrivingVariants();
    }

    @SuppressWarnings("unchecked")
    private void initializeDrivingVariants() {
        // Need to discover the right codec for the driving source of variants manually, since we are
        // treating it specially (separate from the other sources of Features in our FeatureManager).
        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(drivingVariantFile, VariantContext.class);

        //This is the data source for the driving source of variants, which uses a cache lookahead of FEATURE_CACHE_LOOKAHEAD
        drivingVariants = new FeatureDataSource<>(drivingVariantFile, (FeatureCodec<VariantContext, ?>)codec, "drivingVariants", FEATURE_CACHE_LOOKAHEAD);

        //Add the driving datasource to the feature manager too so that it can be queried. Setting lookahead to 0 to avoid caching.
        //Note: we are disabling lookahead here because of windowed queries that need to "look behind" as well.
        drivingVariantsFeatureInput = new FeatureInput<>("drivingVariantFile", Collections.emptyMap(), drivingVariantFile);
        features.addToFeatureSources(0, drivingVariantsFeatureInput, VariantContext.class);

        if ( hasIntervals() ) {
            drivingVariants.setIntervalsForTraversal(intervalsForTraversal);
        }
    }

    /**
     * Returns the feature input for the driving variants file.
     */
    protected final FeatureInput<VariantContext> getDrivingVariantsFeatureInput() {
        return drivingVariantsFeatureInput;
    }

    /**
     * Implementation of variant-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     */
    @Override
    public void traverse() {
        final VariantFilter filter = makeVariantFilter();
        // Process each variant in the input stream.
        StreamSupport.stream(drivingVariants.spliterator(), false)
                .filter(filter)
                .forEach(variant -> {
                    final SimpleInterval variantInterval = new SimpleInterval(variant);
                    apply(variant,
                          new ReadsContext(reads, variantInterval),
                          new ReferenceContext(reference, variantInterval),
                          new FeatureContext(features, variantInterval));

                    progressMeter.update(variantInterval);
                });
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    public final VCFHeader getHeaderForVariants() {
        final Object header = drivingVariants.getHeader();

        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + drivingVariantFile.getAbsolutePath() + " is not in VCF header format");
        }

        return (VCFHeader)header;
    }

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
     *                         an empty array/iterator). Can request extra bases of context around the current read's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext Features spanning the current variant. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     */
    public abstract void apply( VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext );

    /**
     * Close the reads and reference data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalSuccess() instead.
     */
    @Override
    protected final void onShutdown() {
        super.onShutdown();

        if ( drivingVariants != null )
            drivingVariants.close();
    }
}
