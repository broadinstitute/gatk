package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * A Spark version of {@link AssemblyRegionWalker}. Subclasses should implement {@link #processAssemblyRegions(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 */
public abstract class AssemblyRegionWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public final AssemblyRegionReadShardArgumentCollection shardingArgs = new AssemblyRegionReadShardArgumentCollection();

    @ArgumentCollection
    public final AssemblyRegionArgumentCollection assemblyRegionArgs = new AssemblyRegionArgumentCollection();

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    @Argument(doc = "whether to use the strict implementation or not (defaults to the faster implementation that doesn't strictly match the walker version)", fullName = "strict", optional = true)
    public boolean strict = false;

    private String referenceFileName;

    @Override
    public final boolean requiresReads() { return true; }

    @Override
    public final boolean requiresReference() { return true; }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>(2);
        defaultFilters.add(new WellformedReadFilter());
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }

    /**
     * @return The evaluator to be used to determine whether each locus is active or not. Must be implemented by tool authors.
     *         The results of this per-locus evaluator are used to determine the bounds of each active and inactive region.
     */
    public abstract AssemblyRegionEvaluator assemblyRegionEvaluator();

    /**
     * Tools that use an evaluator that is expensive to create, and/or that is not compatible with Spark broadcast, can
     * override this method to return a broadcast of a supplier of the evaluator. The supplier will be invoked once for
     * each Spark partition, thus each partition will have its own evaluator instance.
     */
    protected Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcast(final JavaSparkContext ctx) {
        return assemblyRegionEvaluatorSupplierBroadcastFunction(ctx, assemblyRegionEvaluator());
    }

    private static Broadcast<Supplier<AssemblyRegionEvaluator>> assemblyRegionEvaluatorSupplierBroadcastFunction(final JavaSparkContext ctx, final AssemblyRegionEvaluator assemblyRegionEvaluator) {
        Supplier<AssemblyRegionEvaluator> supplier = (Supplier<AssemblyRegionEvaluator> & Serializable) (() -> assemblyRegionEvaluator);
        return ctx.broadcast(supplier);
    }

    private List<ShardBoundary> intervalShards;

    /**
     * Note that this sets {@code intervalShards} as a side effect, in order to add padding to the intervals.
     */
    @Override
    protected List<SimpleInterval> editIntervals(List<SimpleInterval> rawIntervals) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = rawIntervals == null ? IntervalUtils.getAllIntervalsForReference(sequenceDictionary) : rawIntervals;
        intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardingArgs.readShardSize, shardingArgs.readShardPadding, sequenceDictionary).stream())
                .collect(Collectors.toList());
        List<SimpleInterval> paddedIntervalsForReads =
                intervals.stream().map(interval -> interval.expandWithinContig(shardingArgs.readShardPadding, sequenceDictionary)).collect(Collectors.toList());
        return paddedIntervalsForReads;
    }

    /**
     * Loads assembly regions and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the assembly regions.
     *
     * @return all assembly regions as a {@link JavaRDD}, bounded by intervals if specified.
     */
    protected JavaRDD<AssemblyRegionWalkerContext> getAssemblyRegions(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (strict) {
            return FindAssemblyRegionsSpark.getAssemblyRegionsStrict(ctx, getReads(), getHeaderForReads(), sequenceDictionary, referenceFileName, features,
                    intervalShards, assemblyRegionEvaluatorSupplierBroadcast(ctx), shardingArgs, assemblyRegionArgs,
                    shuffle);
        } else {
            return FindAssemblyRegionsSpark.getAssemblyRegionsFast(ctx, getReads(), getHeaderForReads(), sequenceDictionary, referenceFileName, features,
                    intervalShards, assemblyRegionEvaluatorSupplierBroadcast(ctx), shardingArgs, assemblyRegionArgs,
                    shuffle);
        }
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        referenceFileName = addReferenceFilesForSpark(ctx, referenceArguments.getReferencePath());
        processAssemblyRegions(getAssemblyRegions(ctx), ctx);
    }

    /**
     * Process the assembly regions and write output. Must be implemented by subclasses.
     *
     * @param rdd a distributed collection of {@link AssemblyRegionWalkerContext}
     * @param ctx our Spark context
     */
    protected abstract void processAssemblyRegions(JavaRDD<AssemblyRegionWalkerContext> rdd, JavaSparkContext ctx);

}
