package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.tools.examples.ExampleReferenceWalker;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;

/**
 * A reference walker is a tool which processes each base in a given reference.  Each base can be processed individually
 * or in a context of multiple bases in a window.  Reads and Feature tracks can optionally be included as well. The
 * reference bases to process can be subset to a given set of intervals.
 *
 *  ReferenceWalker authors must implement the apply() method to process each position, and may optionally implement
 *  {@link #onTraversalStart()} and/or {@link #onTraversalSuccess()}.
 *
 *  Tool authors may also implement {@link #getReferenceWindow(SimpleInterval)} to provide additional bases of context
 *  around each position.
 *
 *  See the {@link ExampleReferenceWalker} walker for an example.
 */
public abstract class ReferenceWalker extends WalkerBase {

    @Override
    public String getProgressMeterRecordLabel() { return "bases"; }

    @Override
    public final boolean requiresReference() { return true; }

    /**
     * Initialize data sources for traversal.
     *
     * Marked final so that tool authors don't override it. Tool authors should override onTraversalStart() instead.
     */
    @Override
    protected final void onStartup() {
        super.onStartup();
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of reference-locus-based traversal.
     *
     * NOTE: You should only override {@link #traverse()} if you are writing a new walker base class in the
     * engine package that extends this class. It is not meant to be overridden by tools outside of the engine
     * package.
     */
    @Override
    public void traverse() {
        final CountingReadFilter readFilter = makeReadFilter();

        for(final SimpleInterval locus : getIntervalIterator()){
            final SimpleInterval referenceWindow = getReferenceWindow(locus);
            final ReferenceContext referenceContext = new ReferenceContext(reference, locus, referenceWindow);
            apply(referenceContext,
                  new ReadsContext(reads, referenceContext.getWindow(), readFilter), // Will create an empty ReadsContext if reads == null
                  new FeatureContext(features, referenceContext.getWindow()));   // Will create an empty FeatureContext if features == null

            progressMeter.update(referenceContext.getInterval());
        };

    }

    /**
     * Determine the window to use when creating the ReferenceContext in apply.  This determines which reference bases are
     * seen at each position, as well as which reads / features are considered to overlap the reference site.
     *
     * The default implementation returns the single reference locus passed that is passed in, but subclasses may override
     * this method in order to change the window size.
     *
     * @param locus the current locus being processed in the traversal
     * @return the window to use when creating the ReferenceContext that is passed into {@link #apply(ReferenceContext, ReadsContext, FeatureContext)}
     */
    protected SimpleInterval getReferenceWindow(SimpleInterval locus){
        return locus;
    }

    private Iterable<SimpleInterval> getIntervalIterator(){
        return () -> new IntervalLocusIterator(getTraversalIntervals().iterator());
    }

    /**
     * Process an individual reference locus (with optional contextual information). Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param referenceContext Reference bases at the current locus.  The window around the single locus is provided by {@link #getReferenceWindow(SimpleInterval)}
     *                         Additional bases my be retrieved by invoking the by invoking {@link ReferenceContext#setWindow}
     *                         on this object before calling {@link ReferenceContext#getBases}, however this will not effect
     *                         the set of reads and features that are considered overlapping which are based on the preset
     *                         window passed in.
     * @param readsContext Reads spanning the current reference window. Will be an empty, but non-null, context object
     *                     if there is no backing source of Read data (in which case all queries on it will return an
     *                     empty List).
     * @param featureContext Features spanning the current reference window. Will be an empty, but non-null, context object
     *                       if there is no backing source of Feature data (in which case all queries on it will return an
     *                       empty List).
     *
     */
    public abstract void apply(ReferenceContext referenceContext, ReadsContext readsContext, FeatureContext featureContext );

    /**
     * Shutdown data sources.
     *
     * Marked final so that tool authors don't override it. Tool authors should override {@link #onTraversalSuccess()} instead.
     */
    @Override
    protected final void onShutdown() {
        // Overridden only to make final so that concrete tool implementations don't override
        super.onShutdown();
    }
}
