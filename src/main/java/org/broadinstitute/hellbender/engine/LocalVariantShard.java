package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.iterators.FilteringIterator;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A class to represent a shard of variant data, optionally expanded by a configurable amount of padded data.
 *
 * The reads are lazily loaded by default (when accessing the reads via {@link #iterator}. Loading all the
 * reads in the shard at once is possible via {@link #loadAllRecords()}.
 *
 * The variants returned will overlap the expanded padded interval. It's possible to query whether they are within
 * the main part of the shard via {@link #contains} and {@link #containsStartPosition}.
 *
 * The variants in the shard can be filtered via {@link #loadAllRecords()} (no filtering is performed by default).
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 */
public class LocalVariantShard implements Shard<VariantContext> {

    private final SimpleInterval interval;
    private final SimpleInterval paddedInterval;
    private final GATKDataSource<VariantContext> variantSource;
    private VariantFilter variantFilter;

    /**
     * Create a new Shard spanning the specified interval, with the specified amount of padding.
     *
     * @param interval       the genomic span covered by this shard
     * @param paddedInterval the span covered by this shard, plus any additional padding on each side (must contain the un-padded interval)
     * @param variantSource  source of variants from which to populate this shard
     */
    public LocalVariantShard(final SimpleInterval interval, final SimpleInterval paddedInterval, final GATKDataSource<VariantContext> variantSource) {
        Utils.nonNull(interval);
        Utils.nonNull(paddedInterval);
        Utils.nonNull(variantSource);
        Utils.validateArg(paddedInterval.contains(interval), "The padded interval must contain the un-padded interval");

        this.interval = interval;
        this.paddedInterval = paddedInterval;
        this.variantSource = variantSource;
    }

    /**
     * Create a new Shard spanning the specified interval, with no additional padding
     *
     * @param interval      the genomic span covered by this shard
     * @param variantSource source of reads from which to populate this shard
     */
    public LocalVariantShard(final SimpleInterval interval, final GATKDataSource<VariantContext> variantSource) {
        this(interval, interval, variantSource);
    }

    /**
     * Variants in this shard will be filtered using this filter before being returned.
     *
     * @param filter filter to use (may be null, which signifies that no filtering is to be performed)
     */
    public void setVariantFilter(final VariantFilter filter) {
        this.variantFilter = filter;
    }

    @Override
    public SimpleInterval getInterval() {
        return interval;
    }

    @Override
    public SimpleInterval getPaddedInterval() {
        return paddedInterval;
    }

    /**
     * @return an iterator over variants in this shard, as filtered using the configured variant filter
     * variants are lazily loaded rather than pre-loaded
     */
    @Override
    public Iterator<VariantContext> iterator() {
        final Iterator<VariantContext> iterator = variantSource.query(paddedInterval);
        return (variantFilter == null) ? iterator : new FilteringIterator<VariantContext>(iterator, variantFilter);
    }

    /**
     * Divide an interval into LocalVariantShard. Each shard will cover up to shardSize bases, include shardPadding
     * bases of extra padding on either side, and begin shardStep bases after the previous shard.
     *
     * @param interval interval to shard; must be on the contig according to the provided dictionary
     * @param shardSize desired shard size; intervals larger than this will be divided into shards of up to this size
     * @param shardStep each shard will begin this many bases away from the previous shard
     * @param shardPadding desired shard padding; each shard's interval will be padded on both sides by this number of bases (may be 0)
     * @param variantSource data source for variants
     * @param dictionary sequence dictionary for variants
     * @return List of {@link LocalReadShard} objects spanning the interval
     */
    public static List<LocalVariantShard> divideIntervalIntoShards(final SimpleInterval interval, final int shardSize, final int shardStep, final int shardPadding, final GATKDataSource<VariantContext> variantSource, final SAMSequenceDictionary dictionary) {
        Utils.nonNull(variantSource);
        return Shard.divideIntervalIntoShards(interval, shardSize, shardStep, shardPadding, dictionary)
                .stream().map(shardBoundary -> new LocalVariantShard(shardBoundary.getInterval(), shardBoundary.getPaddedInterval(), variantSource))
                .collect(Collectors.toList());
    }

}
