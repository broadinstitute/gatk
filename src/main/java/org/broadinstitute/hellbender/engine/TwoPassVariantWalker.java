package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.stream.StreamSupport;

/**
 * A modification of a VariantWalker that makes two passes through the variants.
 * This allows the user to store internal states during the first pass, which the user then
 * can process and access during the second pass
 **/
public abstract class TwoPassVariantWalker extends VariantWalker {
    /**
     * Overrides the traversal framework of {@link VariantWalkerBase}
     */
    @Override
    public void traverse(){
        final VariantFilter variantContextFilter = makeVariantFilter();
        final CountingReadFilter readFilter = makeReadFilter();

        // First pass through the variant
        traverseVariants(variantContextFilter, readFilter, this::firstPassApply);
        logger.info("Finished first pass through the variants");

        // Process the data accumulated during the first pass
        afterFirstPass();

        // Second pass
        logger.info("Starting second pass through the variants");
        traverseVariants(variantContextFilter, readFilter, this::secondPassApply);
        logger.info(readFilter.getSummaryLine());
    }

    /**
     *
     * First pass through the variants. Store relevant data as a private member of the walker
     * and process it in {@link #afterFirstPass} before making a second pass
     *
     * @param variant A variant record in vcf
     * @param readsContext Reads overlapping the current variant. Will be empty if a read source (e.g. bam) isn't provided
     * @param referenceContext Reference bases spanning the current variant
     * @param featureContext A vcf record overlapping the current variant from an auxiliary data source (e.g. gnomad)
     */
    protected abstract void firstPassApply(final VariantContext variant,
                                           final ReadsContext readsContext,
                                           final ReferenceContext referenceContext,
                                           final FeatureContext featureContext );

    /**
     *
     * Having seen all of the variants in a vcf, make a second pass through the variants
     *
     * @param variant A variant record in vcf
     * @param readsContext Reads overlapping the current variant. Will be empty if a read source (e.g. bam) isn't provided
     * @param referenceContext Reference bases spanning the current variant
     * @param featureContext A vcf record overlapping the current variant from an auxiliary data source (e.g. gnomad)
     */
    protected abstract void secondPassApply(final VariantContext variant,
                                            final ReadsContext readsContext,
                                            final ReferenceContext referenceContext,
                                            final FeatureContext featureContext);

    /**
     * Process the data collected during the first pass
     */
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

    @FunctionalInterface
    private interface VariantConsumer {
        void consume(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext reference, final FeatureContext features);
    }


    /**
     * Make final to hide it from subclasses
     */
    @Override
    public final void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {}
}
