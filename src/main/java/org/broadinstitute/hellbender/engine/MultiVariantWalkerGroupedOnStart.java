package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


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
    private int firstCurrentVariantStart;
    private int lastCurrentVariantStart;
    private List<ReadsContext> currentReadsContexts = new ArrayList<>();
    private OverlapDetector<SimpleInterval> overlapDetector;

    public static final String IGNORE_VARIANTS_THAT_START_OUTSIDE_INTERVAL = "ignore-variants-starting-outside-interval";

    public static final String COMBINE_VARIANTS_DISTANCE = "combine-variants-distance";

    public static final String MAX_COMBINED_DISTANCE = "max-distance";

    public static final String REFERENCE_WINDOW_PADDING = "ref-padding";

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

    @Advanced
    @Argument(fullName = COMBINE_VARIANTS_DISTANCE, doc = "Maximum distance for variants to be grouped together", optional = true)
    protected int distanceToCombineVariants = defaultDistanceToGroupVariants();

    @Advanced
    @Argument(fullName = MAX_COMBINED_DISTANCE, doc = "Maximum distance for variants to be grouped together", optional = true)
    protected int maxCombinedDistance = defaultMaxGroupedSpan();

    @Advanced
    @Argument(fullName = REFERENCE_WINDOW_PADDING, doc = "Number of bases on either side to expand spanning reference window", optional = true)
    protected int referenceWindowPadding = defaultReferenceWindowPadding();

    // override to group variants that start nearby but not at the same locus
    protected int defaultDistanceToGroupVariants() {
        return 0;
    }

    // override to change reference padding
    protected int defaultReferenceWindowPadding() {
        return 1;
    }

    // override to cap the size span of combined variants
    protected int defaultMaxGroupedSpan() {
        return Integer.MAX_VALUE;
    }

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
            firstCurrentVariantStart = variant.getStart();
        } else if (!currentVariants.get(0).contigsMatch(variant)
                || lastCurrentVariantStart < variant.getStart() - distanceToCombineVariants
                || firstCurrentVariantStart < variant.getStart() - maxCombinedDistance) {
            // Emptying any sites which should emit a new VC since the last one
            apply(new ArrayList<>(currentVariants), currentReadsContexts);
            currentVariants.clear();
            currentReadsContexts.clear();
            firstCurrentVariantStart = variant.getStart();
        }

        currentVariants.add(variant);
        currentReadsContexts.add(readsContext);
        lastCurrentVariantStart = variant.getStart();
    }

    /**
     * This method must be implemented by tool authors.
     *
     * This is the primary traversal for any MultiVariantWalkerGroupedOnStart walkers. Will traverse over input variant contexts
     * and call #apply() exactly once for each unique reference start position. All variants starting at each locus
     * across source files will be grouped and passed as a list of VariantContext objects.
     *  @param variantContexts  VariantContexts from driving variants with matching start position
     *                         NOTE: This will never be empty
     * @param referenceContext  ReferenceContext object covering the reference of the longest spanning VariantContext
     * @param readsContexts
     */
    public abstract void apply(final List<VariantContext> variantContexts, final ReferenceContext referenceContext, final List<ReadsContext> readsContexts);

    public void apply(List<VariantContext> variantContexts, final List<ReadsContext> readsContexts) {
        apply(variantContexts, makeSpanningReferenceContext(variantContexts, referenceWindowPadding), readsContexts);
    }

    /**
     * Helper method that ensures the reference context it returns is adequate to span the length of all the accumulated
     * VariantContexts. It assumes that all variant contexts in currentVariants have the same contig.
     */
    private ReferenceContext makeSpanningReferenceContext(final List<VariantContext> variantContexts, final int referenceWindowPadding) {
        Utils.nonEmpty(variantContexts, "Must have at least one current variant context");
        final List<String> contigs = variantContexts.stream().map(VariantContext::getContig).distinct().collect(Collectors.toList());
        Utils.validate(contigs.size() == 1, "variant contexts should all have the same contig");
        final int minStart = variantContexts.stream().mapToInt(VariantContext::getStart).min().getAsInt();
        final int maxEnd = variantContexts.stream().mapToInt(VariantContext::getEnd).max().getAsInt();
        final SimpleInterval combinedInterval = new SimpleInterval(contigs.get(0), minStart, maxEnd);

        final ReferenceContext combinedReferenceContext = new ReferenceContext(reference, combinedInterval);
        combinedReferenceContext.setWindow(referenceWindowPadding,referenceWindowPadding);
        return combinedReferenceContext;
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of multi-variant grouped on start traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the
     * engine package.
     */
    @Override
    public void traverse() {
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
            apply(currentVariants, currentReadsContexts);
        }
    }
}
