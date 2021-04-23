package org.broadinstitute.hellbender.engine.spark;

import com.google.common.base.Function;
import com.google.common.collect.Iterators;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import javax.annotation.Nullable;
import java.util.Iterator;

/**
 * A Spark version of {@link ReadWalker}. Subclasses should implement {@link #processReads(JavaRDD, JavaSparkContext)}
 * and operate on the passed in RDD.
 * @see org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithReferenceSpark
 * @see org.broadinstitute.hellbender.tools.examples.ExampleReadWalkerWithVariantsSpark
 */
public abstract class ReadWalkerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() {
        return true;
    }

    /**
     * This number controls the size of the cache for our FeatureInputs
     * (specifically, the number of additional bases worth of overlapping records to cache when querying feature sources).
     */
    public static final int FEATURE_CACHE_LOOKAHEAD = 1_000;

    private String referenceFileName;

    void initializeFeatures() {
        features = new FeatureManager(this, FEATURE_CACHE_LOOKAHEAD);
        if ( features.isEmpty() ) {  // No available sources of Features discovered for this tool
            features = null;
        }
    }

    /**
     * Loads reads and the corresponding reference and features into a {@link JavaRDD} for the intervals specified.
     *
     * If no intervals were specified, returns all the reads.
     *
     * @return all reads as a {@link JavaRDD}, bounded by intervals if specified.
     */
    public JavaRDD<ReadWalkerContext> getReads(JavaSparkContext ctx) {
        Broadcast<FeatureManager> bFeatureManager = features == null ? null : ctx.broadcast(features);
        return getReads().mapPartitions(getReadsFunction(referenceFileName, bFeatureManager));
    }

    private static FlatMapFunction<Iterator<GATKRead>, ReadWalkerContext> getReadsFunction(
            String referenceFileName, Broadcast<FeatureManager> bFeatureManager) {
        return readIterator -> {
            ReferenceDataSource reference = referenceFileName == null ? null : new ReferenceFileSource(IOUtils.getPath(SparkFiles.get(referenceFileName)));
            FeatureManager features = bFeatureManager == null ? null : bFeatureManager.getValue();
            return Iterators.transform(readIterator, new Function<GATKRead, ReadWalkerContext>() {
                @Nullable
                @Override
                public ReadWalkerContext apply(@Nullable GATKRead r) {
                    final SimpleInterval readInterval = getReadInterval(r);
                    return new ReadWalkerContext(r, new ReferenceContext(reference, readInterval), new FeatureContext(features, readInterval));
                }
            });
        };
    }

    /**
     * Returns an interval for the read.
     * Note: some walkers must be able to work on any read, including those whose coordinates do not form a valid SimpleInterval.
     * So here we check this condition and create null intervals for such reads.
     */
    static SimpleInterval getReadInterval(final GATKRead read) {
        return !read.isUnmapped() && SimpleInterval.isValid(read.getContig(), read.getStart(), read.getEnd()) ? new SimpleInterval(read) : null;
    }

    @Override
    protected void runTool(JavaSparkContext ctx) {
        referenceFileName = addReferenceFilesForSpark(ctx, referenceArguments.getReferencePath());
        processReads(getReads(ctx), ctx);
    }

    /**
     * Process the reads and write output. Must be implemented by subclasses.
     *
     * @param rdd a distributed collection of {@link ReadWalkerContext}
     * @param ctx our Spark context
     */
    protected abstract void processReads(JavaRDD<ReadWalkerContext> rdd, JavaSparkContext ctx);
}
