package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ListMap;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.*;

/**
 * Iterator that dynamically applies filter strings to VariantContext records supplied by an underlying
 * iterator.  Returns all records from the underlying stream and does not remove any.
 *
 * @author tfennell
 */
public final class FilterApplyingVariantIterator implements CloseableIterator<VariantContext> {
    /** Filter string that is used to filter a Variant when all variant genotypes are filtered out. */
    public static final String ALL_GTS_FILTERED = "AllGtsFiltered";
    /** The "PASS"ing filter String. */
    public static final String PASS_FILTER = "PASS";

    private final Iterator<VariantContext> iterator;
    private final VariantFilter[] filters;
    private final GenotypeFilter[] gtFilters;

    /**
     * Constructs an iterator from an underlying iterator and the provided (possibly empty)
     * collections of variant and genotype filters.
     */
    public FilterApplyingVariantIterator(final Iterator<VariantContext> iterator,
                                         final Collection<VariantFilter> filters,
                                         final Collection<GenotypeFilter> gtFilters) {
        this.iterator = iterator;
        this.filters = filters.toArray(new VariantFilter[filters.size()]);
        this.gtFilters = gtFilters.toArray(new GenotypeFilter[gtFilters.size()]);
    }

    /**
     * Provides the next record from the underlying iterator after applying filter strings generated
     * by the set of filters in use by the iterator.
     */
    @Override
    public VariantContext next() {
        final VariantContext ctx = this.iterator.next();
        final Set<String> filterStrings = new HashSet<>();

        // Collect variant level filters
        for (final VariantFilter filter : this.filters) {
            final String val = filter.filter(ctx);
            if (val != null) filterStrings.add(val);
        }

        // Collect genotype level filters in a Map of Sample -> List<filter string>
        final ListMap<String,String> gtFilterStrings = new ListMap<>();
        final Set<String> variantSamples = new HashSet<>();
        for (final Genotype gt : ctx.getGenotypes()) {
            if (gt.isCalled() && !gt.isHomRef()) variantSamples.add(gt.getSampleName());

            for (final GenotypeFilter filter : gtFilters) {
                final String filterString = filter.filter(ctx,gt);
                if (filterString != null)  gtFilterStrings.add(gt.getSampleName(), filterString);
            }
        }

        // If all genotypes are filtered apply a site level filter
        if (gtFilterStrings.keySet().containsAll(variantSamples)) {
            filterStrings.add(ALL_GTS_FILTERED);
        }

        // Make a builder and set the site level filter appropriately
        final VariantContextBuilder builder = new VariantContextBuilder(ctx);
        if (filterStrings.isEmpty()) {
            builder.passFilters();
        }
        else {
            builder.filters(filterStrings);
        }

        // Apply filters to the necessary genotypes
        builder.noGenotypes();
        final List<Genotype> newGenotypes = new ArrayList<>(ctx.getNSamples());
        for (final Genotype gt : ctx.getGenotypes()) {
            final GenotypeBuilder gtBuilder = new GenotypeBuilder(gt);
            final List<String> filters = gtFilterStrings.get(gt.getSampleName());

            if (filters == null || filters.isEmpty()) {
                gtBuilder.filter(PASS_FILTER);
            }
            else {
                gtBuilder.filters(filters);
            }
            newGenotypes.add(gtBuilder.make());
        }
        builder.genotypes(newGenotypes);

        return builder.make();
    }

    @Override
    public boolean hasNext() { return this.iterator.hasNext(); }
    @Override
    public void close() { CloserUtil.close(this.iterator); }
    @Override
    public void remove() { throw new UnsupportedOperationException("remove() not supported by FilterApplyingVariantIterator."); }
}
