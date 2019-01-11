package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.List;


/**
 * A MultiVariantWalker that walks over multiple variant context sources in reference order and emits to client tools
 * groups of all input variant contexts by their start position. This is intended to mimic GATK3 traversal behavior for
 * some tools.
 *
 * As such, the argument '-ignore-variants-starting-outside-interval' has been provided to mimic GATK3's behavior
 * only presenting variants that start inside the requested interval regardless of whether there is a spanning variant.
 *
 * Client tools must implement apply(List<VariantContext> variantContexts, ReferenceContext referenceContext)
 */
public abstract class MultiVariantWalkerGroupedOnStart extends MultiVariantWalker {
    private List<VariantContext> currentVariants = new ArrayList<>();
    private ReferenceContext spanningReferenceContext;
    private OverlapDetector<SimpleInterval> overlapDetector;

    public static final String IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL = "ignore-variants-starting-outside-interval";

    /**
     * this option has no effect unless intervals are specified.
     * <p>
     * This exists to mimic GATK3 interval traversal patterns
     */
    @Advanced
    @Argument(fullName = IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL,
            doc = "Restrict variant output to sites that start within provided intervals (only applies when an interval is specified)",
            optional = true)
    private boolean ignoreIntervalsOutsideStart = false;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /**
     * This method keeps track of all the variants it is passed and will feed all the variants that start at the same
     * site to the reduce method.
     *
     * {@inheritDoc}
     */
    @Override
    public final void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {

        // Filtering out variants that start outside of the specified intervals
        if (ignoreIntervalsOutsideStart && !isWithinInterval(new SimpleInterval(variant.getContig(), variant.getStart(), variant.getStart()))) {
            return;
        }

        // Collecting all the reads that start at a particular base into one.
        if (currentVariants.isEmpty()) {
            currentVariants.add(variant);
        } else if (!currentVariants.get(0).contigsMatch(variant)
                || currentVariants.get(0).getStart() < variant.getStart()) {
            // Emptying any sites which should emit a new VC since the last one
            apply(new ArrayList<>(currentVariants), spanningReferenceContext);
            currentVariants.clear();
            currentVariants.add(variant);
        } else {
            currentVariants.add(variant);
        }
        if (referenceContext.hasBackingDataSource()){
            referenceContext.setWindow(1, 1);
        }
        spanningReferenceContext = getExpandedReferenceContext(currentVariants, spanningReferenceContext, referenceContext);
    }

    /**
     * This method must be implemented by tool authors.
     *
     * This is the primary traversal for any MultiVariantWalkerGroupedOnStart walkers. Will traverse over input variant contexts
     * and call #apply() exactly once for each unique reference start position. All variants starting at each locus
     * across source files will be grouped and passed as a list of VariantContext objects.
     *
     * @param variantContexts  VariantContexts from driving variants with matching start positon
     *                         NOTE: This will never be empty
     * @param referenceContext  ReferenceContext object covering the reference of the longest spanning VariantContext
     *                          with one base on either end as window padding.
     */
    public abstract void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext);

    /**
     * Helper method that ensures the reference context it returns is adequate to span the length of all the accumulated
     * VariantContexts. It makes the assumption that all variant contexts in currentVariants start at the same location and
     * that currentReferenceContext corresponds to the correct span for the longest variantContexts the first N-1 VariantContexts
     * and that newReferenceContext corresponds to the correct span for the Nth item.
     */
    @VisibleForTesting
    final static ReferenceContext getExpandedReferenceContext(List<VariantContext> currentVariants,
                                                              ReferenceContext currentReferenceContext,
                                                              ReferenceContext newReferenceContext) {
        if ((currentReferenceContext==null)||(currentVariants.size() == 1 || newReferenceContext.getWindow().getEnd() > currentReferenceContext.getWindow().getEnd())) {
            return newReferenceContext;
        }
        return currentReferenceContext;
    }

    @Override
    public final void traverse() {
        beforeTraverse();
        super.traverse();
        afterTraverse();
    }

    /**
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalStart} instead.
     */
    private void beforeTraverse() {
        overlapDetector = hasUserSuppliedIntervals() ? OverlapDetector.create(intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary())) : null;
    }

    /**
     * @param loc locatable to query
     * @return true if the query loc is entirely contained by the interval, true if no interval
     */
    protected final boolean isWithinInterval(Locatable loc) {
        return (overlapDetector==null || overlapDetector.overlapsAny(loc));
    }

    /**
     * Clear accumulated reads before {@link #onTraversalSuccess()} is accessed
     */
    private void afterTraverse() {
        // Clearing the accumulator
        if (currentVariants.isEmpty()) {
            logger.warn("Error: The requested interval contained no data in source VCF files");

        } else {
            apply(currentVariants, spanningReferenceContext);
        }
    }
}
