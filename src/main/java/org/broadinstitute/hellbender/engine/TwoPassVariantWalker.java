package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.stream.StreamSupport;

public abstract class TwoPassVariantWalker extends VariantWalker {

    @Override
    public void traverse(){
        final VariantFilter variantFilter = makeVariantFilter();
        final CountingReadFilter readFilter = makeReadFilter();

        // First pass through the variant
        traverseVariants(variantFilter, readFilter, this::firstPassApply);
        logger.info("Finished first pass through the variants");

        afterFirstPass();

        // Second pass
        logger.info("Starting second pass through the variants");
        traverseVariants(variantFilter, readFilter, this::secondPassApply);
        logger.info(readFilter.getSummaryLine());
    }


    protected abstract void firstPassApply(final VariantContext variant,
                                           final ReadsContext readsContext,
                                           final ReferenceContext referenceContext,
                                           final FeatureContext featureContext );

    protected abstract void secondPassApply(final VariantContext variant,
                                            final ReadsContext readsContext,
                                            final ReferenceContext referenceContext,
                                            final FeatureContext featureContext);

    protected abstract void afterFirstPass();

    private void traverseVariants(final VariantFilter variantFilter, final CountingReadFilter readFilter, final VariantConsumer variantConsumer){
        StreamSupport.stream(getSpliteratorForDrivingVariants(), false)
                .filter(variantFilter)
                .forEach(variant -> {
                    final SimpleInterval variantInterval = new SimpleInterval(variant);
                    variantConsumer.consume(variant,
                            new ReadsContext(reads, variantInterval, readFilter),
                            new ReferenceContext(reference, variantInterval),
                            new FeatureContext(features, variantInterval));

                    progressMeter.update(variantInterval);
                });
    }

    /**
     * a common abstraction for first and second pass apply functions
     */
    @FunctionalInterface
    private interface VariantConsumer {
        void consume(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext reference, final FeatureContext features);
    }


    // Make this final to hide it from subclasses
    @Override
    public final void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
    }

}
