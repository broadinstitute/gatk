package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Spliterator;

/**
 * A VariantWalker is a tool that processes a variant at a time from a source of variants, with
 * optional contextual information from a reference, sets of reads, and/or supplementary sources
 * of Features.
 *
 * VariantWalker authors must implement the {@link #apply} method to process each variant, and may optionally implement
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class VariantWalker extends VariantWalkerBase {

    // NOTE: using GATKPath rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public GATKPath drivingVariantFile;

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager in GATKTool
    // we do add the driving source to the Feature manager but we do need to treat it differently and thus this field.
    private FeatureDataSource<VariantContext> drivingVariants;
    private FeatureInput<VariantContext> drivingVariantsFeatureInput;

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

    @Override
    protected void initializeDrivingVariants() {
        drivingVariantsFeatureInput = new FeatureInput<>(drivingVariantFile, "drivingVariantFile");

        // Create a FeatureDataSource for the driving variants FeatureInput, using the
        // cache lookahead value from getDrivingVariantCacheLookAheadBases()
        drivingVariants = new FeatureDataSource<>(drivingVariantsFeatureInput, getDrivingVariantCacheLookAheadBases(), VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                                  getGenomicsDBOptions());

        // Also add the driving variants FeatureInput to FeatureManager as well so that it can be queried,
        // but use a lookahead value of 0 to avoid caching because of windowed queries that need to "look behind" as well.
        features.addToFeatureSources(0, drivingVariantsFeatureInput, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                                     getGenomicsDBOptions());

        // Note: the intervals for the driving variants are set in onStartup()
    }

    /**
     * Returns the feature input for the driving variants file.
     */
    protected final FeatureInput<VariantContext> getDrivingVariantsFeatureInput() {
        return drivingVariantsFeatureInput;
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    public final VCFHeader getHeaderForVariants() {
        final Object header = drivingVariants.getHeader();

        if ( ! (header instanceof VCFHeader) ) {
            throw new GATKException("Header for " + drivingVariantFile + " is not in VCF header format");
        }

        return (VCFHeader)header;
    }

    /**
     * Implementation of variant-based traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the engine
     * package.
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

        if ( drivingVariants != null )
            drivingVariants.close();
    }
}
