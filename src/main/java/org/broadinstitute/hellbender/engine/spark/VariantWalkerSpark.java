package org.broadinstitute.hellbender.engine.spark;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.VariantFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.utils.IndexUtils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple4;

import javax.annotation.Nullable;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A Spark version of {@link VariantWalker}.
 */
public abstract class VariantWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    // NOTE: using File rather than FeatureInput<VariantContext> here so that we can keep this driving source
    //       of variants separate from any other potential sources of Features
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
    public String drivingVariantFile;

    @Argument(fullName="variantShardSize", shortName="variantShardSize", doc = "Maximum size of each variant shard, in bases.", optional = true)
    public int variantShardSize = 10000;

    @Argument(fullName="variantShardPadding", shortName="readShardPadding", doc = "Each variant shard has this many bases of extra context on each side.", optional = true)
    public int variantShardPadding = 1000;

    @Argument(doc = "whether to use the shuffle implementation or not", shortName = "shuffle", fullName = "shuffle", optional = true)
    public boolean shuffle = false;

    private transient VariantsSparkSource variantsSource;

    /**
     * This number controls the size of the cache for our primary and auxiliary FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 100_000;

    @Override
    protected void runPipeline(JavaSparkContext sparkContext) {
        initializeVariants(sparkContext);
        super.runPipeline(sparkContext);
    }

    void initializeVariants(final JavaSparkContext sparkContext) {
        variantsSource = new VariantsSparkSource(sparkContext);
    }

    void initializeFeatures() {
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    @Override
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {
        final SAMSequenceDictionary dictFromDrivingVariants = getHeaderForVariants().getSequenceDictionary();
        if (dictFromDrivingVariants != null){
            //If this dictionary looks like it was synthesized from a feature index, see
            //if there is a better dictionary available from another source (i.e., the reference)
            if (IndexUtils.isSequenceDictionaryFromIndex(dictFromDrivingVariants)) {
                final SAMSequenceDictionary otherDictionary = super.getBestAvailableSequenceDictionary();
                return otherDictionary != null ?
                        otherDictionary :
                        dictFromDrivingVariants;
            } else {
                return dictFromDrivingVariants;
            }
        }
        return super.getBestAvailableSequenceDictionary();
    }

    /**
     * Gets the header associated with our driving source of variants as a VCFHeader.
     *
     * @return VCFHeader for our driving source of variants
     */
    public final VCFHeader getHeaderForVariants() {
        return VariantsSparkSource.getHeader(drivingVariantFile);
    }

    /**
     * Returns the variant filter (simple or composite) that will be applied to the variants before calling {@link #getVariants}.
     * The default implementation filters nothing.
     * Default implementation of {@link #getVariants} calls this method once before iterating
     * over the reads and reuses the filter object to avoid object allocation. Nevertheless, keeping state in filter objects is strongly discouraged.
     *
     * Subclasses can extend to provide own filters (ie override and call super).
     * Multiple filters can be composed by using {@link VariantFilter} composition methods.
     */
    protected VariantFilter makeVariantFilter() {
        return VariantFilterLibrary.ALLOW_ALL_VARIANTS;
    }

    /**
     * Loads variants and the corresponding reads, reference and features into a {@link JavaRDD} for the intervals specified.
     * FOr the current implementation the reads context will always be empty.
     *
     * If no intervals were specified, returns all the variants.
     *
     * @return all variants as a {@link JavaRDD}, bounded by intervals if specified.
     */
    public JavaRDD<Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext>> getVariants(JavaSparkContext ctx) {
        SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> intervalShards = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, variantShardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        JavaRDD<VariantContext> variants = variantsSource.getParallelVariantContexts(drivingVariantFile, getIntervals());
        VariantFilter variantFilter = makeVariantFilter();
        variants = variants.filter(variantFilter::test);
        JavaRDD<Shard<VariantContext>> shardedVariants = SparkSharder.shard(ctx, variants, VariantContext.class, sequenceDictionary, intervalShards, variantShardSize, shuffle);
        Broadcast<ReferenceMultiSource> bReferenceSource = hasReference() ? ctx.broadcast(getReference()) : null;
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return shardedVariants.flatMap(getVariantsFunction(bReferenceSource, bFeatureManager, sequenceDictionary, variantShardPadding));
    }

    private static FlatMapFunction<Shard<VariantContext>, Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext>> getVariantsFunction(
            final Broadcast<ReferenceMultiSource> bReferenceSource,
            final Broadcast<FeatureManager> bFeatureManager,
            final SAMSequenceDictionary sequenceDictionary, final int variantShardPadding) {
        return (FlatMapFunction<Shard<VariantContext>, Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext>>) shard -> {
            // get reference bases for this shard (padded)
            SimpleInterval paddedInterval = shard.getInterval().expandWithinContig(variantShardPadding, sequenceDictionary);
            ReferenceDataSource reference = bReferenceSource == null ? null :
                    new ReferenceMemorySource(bReferenceSource.getValue().getReferenceBases(null, paddedInterval), sequenceDictionary);
            FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();

            Iterator<Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext>> transform = Iterators.transform(shard.iterator(), new Function<VariantContext, Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext>>() {
                @Nullable
                @Override
                public Tuple4<VariantContext, ReadsContext, ReferenceContext, FeatureContext> apply(@Nullable VariantContext variant) {
                    final SimpleInterval variantInterval = new SimpleInterval(variant);
                    return new Tuple4<>(variant,
                            new ReadsContext(), // empty
                            new ReferenceContext(reference, variantInterval),
                            new FeatureContext(features, variantInterval));
                }
            });
            // only include variants that start in the shard
            return Iterators.filter(transform, r -> r._1().getStart() >= shard.getStart()
                    && r._1().getStart() <= shard.getEnd());
        };
    }
}
