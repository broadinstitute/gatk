package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.transformers.VariantTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * VariantWalkerGroupedByLocus processes variants from a single source, either one at a time in order, or optionally
 * grouped by locus, with optional contextual information from a reference, sets of reads, and/or supplementary sources
 * of Features. By-locus traversal is opt in, via {@link #setGroupByLocus()}, and only loci with overlapping
 * variants are traversed.
 *
 * VariantWalkerGroupedByLocus authors must implement the {@link #apply} method to process each variant, and may optionally implement
 * {@link #onTraversalStart}, {@link #onTraversalSuccess} and/or {@link #closeTool}.
 */
public abstract class VariantWalkerGroupedByLocus extends VariantWalkerBase {

    // NOTE: using File rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public String drivingVariantFile;

    // NOTE: keeping the driving source of variants separate from other, supplementary FeatureInputs in our FeatureManager in GATKTool
    //we do add the driving source to the Feature manager but we do need to treat it differently and thus this field.
    private FeatureDataSource<VariantContext> drivingVariants;
    private FeatureInput<VariantContext> drivingVariantsFeatureInput;

    private boolean traverseByLocus = false;

    /**
     * Enable grouped of variants by locus. When true, the {@link #apply} method will be called with all overlapping
     * variants at each each locus that has overlapping variants.
     */
    protected void setGroupByLocus() { traverseByLocus = true; }

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

        //This is the data source for the driving source of variants, which uses a cache lookahead of FEATURE_CACHE_LOOKAHEAD
        drivingVariants = new FeatureDataSource<>(drivingVariantsFeatureInput, FEATURE_CACHE_LOOKAHEAD, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                referenceArguments.getReferencePath());

        //Add the driving datasource to the feature manager too so that it can be queried. Setting lookahead to 0 to avoid caching.
        //Note: we are disabling lookahead here because of windowed queries that need to "look behind" as well.
        features.addToFeatureSources(0, drivingVariantsFeatureInput, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer,
                referenceArguments.getReferencePath());

        //Note: the intervals for the driving variants are set in onStartup
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

    @Override
    public boolean requiresReference() {
        return true;
    }

    /**
     * Implementation of variant-based traversal.
     * Subclasses can override to provide their own behavior but default implementation should be suitable for most uses.
     */
    @Override
    public void traverse() {
        final CountingReadFilter readFilter = makeReadFilter();
        final VariantFilter variantFilter = makeVariantFilter();

        if (traverseByLocus) {
            getLocusStream(getTraversalIntervals())
                .forEach(locus -> {
                    // Get all of the variants overlapping this locus, filter them and collect into a list
                    final Iterator<VariantContext> overlappingVariants = drivingVariants.query(locus);
                    final List<VariantContext> filteredVariants = getTransformedVariantStream(
                            Spliterators.spliteratorUnknownSize(overlappingVariants, 0), variantFilter)
                            .collect(Collectors.toList());
                        if (!filteredVariants.isEmpty()) {
                            apply(locus,
                                    filteredVariants,
                                    new ReadsContext(reads, locus, readFilter),
                                    new ReferenceContext(reference, locus),
                                    new FeatureContext(features, locus));

                            progressMeter.update(locus);
                        }
                    }
            );
        } else {
            // Process each variant in the input stream.
            getTransformedVariantStream( getSpliteratorForDrivingVariants(), variantFilter )
                    .forEach(variant -> {
                        final SimpleInterval variantInterval = new SimpleInterval(variant);
                        apply(variant,
                                Collections.singletonList(variant),
                                new ReadsContext(reads, variantInterval, readFilter),
                                new ReferenceContext(reference, variantInterval),
                                new FeatureContext(features, variantInterval));

                        progressMeter.update(variantInterval);
                    });
        }
    }

    protected Stream<VariantContext> getTransformedVariantStream(final Spliterator<VariantContext> source, final VariantFilter filter) {
        final VariantTransformer preTransformer  = makePreVariantFilterTransformer();
        final VariantTransformer postTransformer = makePostVariantFilterTransformer();
        return StreamSupport.stream(source, false)
                .map(preTransformer)
                .filter(filter)
                .map(postTransformer);
    }

    // Return a Stream of SimpleInterval covering the entire territory sketched out by requestedIntervals
    private Stream<SimpleInterval> getLocusStream(final List<SimpleInterval> requestedIntervals) {
        final Iterable<SimpleInterval> iterable = () -> new IntervalLocusIterator(requestedIntervals.iterator());
        return StreamSupport.stream(iterable.spliterator(), false);
    }

    /**
     * Process an individual variant. Must be implemented by tool authors.
     * In general, tool authors should simply stream their output from apply(), and maintain as little internal state
     * as possible.
     *
     * @param loc
     * @param variants Current variants being processed.
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
    public abstract void apply(Locatable loc, List<VariantContext> variants, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext );

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
