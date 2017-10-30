package org.broadinstitute.hellbender.engine;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;


/**
 * A MultiVariantWalker that walks over multiple variant context sources in reference order and emits to client tools
 * groups of all input variant contexts by their start position. This is intended to mimic GATK3 traversal behavior for
 * some tools.
 *
 * As such, the argument '-ignore_variants_starting_outside_interval' has been provided to mimic GATK3's behavior
 * only presenting variant that start inside the requested interval regardless of whether there is a spanning variant.
 *
 * Client tools must implement apply(List<VariantContext> variantContexts, ReferenceContext referenceContext)
 */
public abstract class MultiVariantWalkerGroupedOnStart extends MultiVariantWalker {
    private List<VariantContext> currentVariants = new ArrayList<>();
    private ReferenceContext spanningReferenceContext;
    private OverlapDetector<SimpleInterval> overlapDetector;

    public static final String IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL = "ignore_variants_starting_outside_interval";

    /**
     * This option can only be activated if intervals are specified.
     * <p>
     * This exists to mimic GATK3 interval traversal patterns
     */
    @Advanced
    @Argument(fullName = IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL,
            doc = "Restrict variant output to sites that start within provided intervals",
            optional = true)
    private boolean ignoreIntervalsOutsideStart = false;

    /**
     * This method keeps track of all the variants it is passed and will feed all the variants that start at the same
     * site to the reduce method.
     *
     * @param variant          Current variant being processed.
     * @param readsContext     Reads overlapping the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of reads data (in which case all queries on it will return
     *                         an empty array/iterator)
     * @param referenceContext Reference bases spanning the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of reference data (in which case all queries on it will return
     *                         an empty array/iterator). Can request extra bases of context around the current variant's interval
     *                         by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}
     * @param featureContext   Features spanning the current variant. Will be an empty, but non-null, context object
     *                         if there is no backing source of Feature data (in which case all queries on it will return an
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
            apply(currentVariants, spanningReferenceContext);
            currentVariants.clear();
            currentVariants.add(variant);
        } else {
            currentVariants.add(variant);
        }
        referenceContext.setWindow(1, 1);
        spanningReferenceContext = updatePositionalState(currentVariants, spanningReferenceContext, referenceContext);
    }

    /**
     * Primary traversal, feeds a list of variant contexts across driving variant files with matching start locations in
     * sorted order based on start location.
     *
     * @param variantContexts  VariantContexts from driving variants with matching start positon
     * @param referenceContext  ReferenceContext object covering the reference of the longest spanning VariantContext
     *                          with one base on either end as window padding.
     */
    public abstract void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext);

    /**
     * Method which updates that the currentQueuedContextState object to add a new reference context making sure
     * to hold the longest set of reference bases it has seen, which should all originate from the same source
     *
     * @param currentVariants
     */
    @VisibleForTesting
    final static ReferenceContext updatePositionalState(List<VariantContext> currentVariants,
                                                        ReferenceContext currentReferenceContext,
                                                        ReferenceContext newReferenceContext) {
        if (currentVariants.size() == 1 || newReferenceContext.getWindow().getEnd() > currentReferenceContext.getWindow().getEnd()) {
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
    protected final void beforeTraverse() {
        overlapDetector = hasIntervals() ? OverlapDetector.create(intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary())) : null;
    }

    /**
     * @param loc locatable to query
     * @return true if the query loc is entirely contained by the interval, true if no interval
     */
    public final boolean isWithinInterval(Locatable loc) {
        return (overlapDetector==null || overlapDetector.overlapsAny(loc));
    }

    /**
     * Clear accumulated reads before {@link #onTraversalSuccess()} is accessed
     */
    public void afterTraverse() {
        // Clearing the accumulator
        if (currentVariants.isEmpty()) {
            logger.warn("Error: The requested interval contained no data in source VCF files");

        } else {
            apply(currentVariants, spanningReferenceContext);
        }
    }
}
