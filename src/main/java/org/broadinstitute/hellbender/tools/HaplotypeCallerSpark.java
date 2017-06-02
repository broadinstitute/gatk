package org.broadinstitute.hellbender.tools;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.SparkReadShard;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Call germline SNPs and indels via local re-assembly of haplotypes
 *
 * This is an implementation of {@link HaplotypeCaller} using spark to distribute the computation.
 * It is still in an early stage of development and does not yet support all the options that the non-spark version does.
 *
 * Specifically it does not support the --dbsnp, --comp, and --bamOutput options.
 *
 */
@CommandLineProgramProperties(summary = "HaplotypeCaller on Spark", oneLineSummary = "HaplotypeCaller on Spark", programGroup = SparkProgramGroup.class)
public final class HaplotypeCallerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Single file to which variants should be written")
    public String output;

    @ArgumentCollection
    public final ShardingArgumentCollection shardingArgs = new ShardingArgumentCollection();

    public static class ShardingArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        @Argument(fullName="readShardSize", shortName="readShardSize", doc = "Maximum size of each read shard, in bases. For good performance, this should be much larger than the maximum assembly region size.", optional = true)
        public int readShardSize = HaplotypeCaller.DEFAULT_READSHARD_SIZE;

        @Argument(fullName="readShardPadding", shortName="readShardPadding", doc = "Each read shard has this many bases of extra context on each side. Read shards must have as much or more padding than assembly regions.", optional = true)
        public int readShardPadding = HaplotypeCaller.DEFAULT_READSHARD_PADDING;

        @Argument(fullName = "minAssemblyRegionSize", shortName = "minAssemblyRegionSize", doc = "Minimum size of an assembly region", optional = true)
        public int minAssemblyRegionSize = HaplotypeCaller.DEFAULT_MIN_ASSEMBLY_REGION_SIZE;

        @Argument(fullName = "maxAssemblyRegionSize", shortName = "maxAssemblyRegionSize", doc = "Maximum size of an assembly region", optional = true)
        public int maxAssemblyRegionSize = HaplotypeCaller.DEFAULT_MAX_ASSEMBLY_REGION_SIZE;

        @Argument(fullName = "assemblyRegionPadding", shortName = "assemblyRegionPadding", doc = "Number of additional bases of context to include around each assembly region", optional = true)
        public int  assemblyRegionPadding = HaplotypeCaller.DEFAULT_ASSEMBLY_REGION_PADDING;

        @Advanced
        @Argument(fullName = "activeProbabilityThreshold", shortName = "activeProbabilityThreshold", doc="Minimum probability for a locus to be considered active.", optional = true)
        public double activeProbThreshold = HaplotypeCaller.DEFAULT_ACTIVE_PROB_THRESHOLD;

        @Advanced
        @Argument(fullName = "maxProbPropagationDistance", shortName = "maxProbPropagationDistance", doc="Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions", optional = true)
        public int maxProbPropagationDistance = HaplotypeCaller.DEFAULT_MAX_PROB_PROPAGATION_DISTANCE;

    }

    @ArgumentCollection
    public HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    @Override
    public boolean requiresReads(){
        return true;
    }

    @Override
    public boolean requiresReference(){
        return true;
    }


    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(getHeaderForReads().getSequenceDictionary());
        final JavaRDD<VariantContext> variants = callVariantsWithHaplotypeCaller(getAuthHolder(), ctx, getReads(), getHeaderForReads(), getReference(), intervals, hcArgs, shardingArgs);
        writeVariants(variants);
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return HaplotypeCallerEngine.makeStandardHCReadFilters();
    }

    /**
     * Call Variants using HaplotypeCaller on Spark and return an RDD of  {@link VariantContext}
     *
     * This may be called from any spark pipeline in order to call variants from an RDD of GATKRead
     *
     * @param authHolder authorization needed for the reading the reference
     * @param ctx the spark context
     * @param reads the reads variants should be called from
     * @param header the header that goes with the reads
     * @param reference the reference to use when calling
     * @param intervals the intervals to restrict calling to
     * @param hcArgs haplotype caller arguments
     * @param shardingArgs arguments to control how the assembly regions are sharded
     * @return an RDD of Variants
     */
    public static JavaRDD<VariantContext> callVariantsWithHaplotypeCaller(
            final AuthHolder authHolder,
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads,
            final SAMFileHeader header,
            final ReferenceMultiSource reference,
            final List<SimpleInterval> intervals,
            final HaplotypeCallerArgumentCollection hcArgs,
            final ShardingArgumentCollection shardingArgs) {
        Utils.validateArg(hcArgs.dbsnp.dbsnp == null, "HaplotypeCallerSpark does not yet support -D or --dbsnp arguments" );
        Utils.validateArg(hcArgs.comps.isEmpty(), "HaplotypeCallerSpark does not yet support -comp or --comp arguments" );
        Utils.validateArg(hcArgs.bamOutputPath == null, "HaplotypeCallerSpark does not yet support -bamout or --bamOutput");
        if ( !reference.isCompatibleWithSparkBroadcast()){
            throw new UserException.Require2BitReferenceForBroadcast();
        }

        final Broadcast<ReferenceMultiSource> referenceBroadcast = ctx.broadcast(reference);
        final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast = ctx.broadcast(hcArgs);
        final OverlapDetector<ShardBoundary> overlaps = getShardBoundaryOverlapDetector(header, intervals, shardingArgs.readShardSize, shardingArgs.readShardPadding);
        final Broadcast<OverlapDetector<ShardBoundary>> shardBoundariesBroadcast = ctx.broadcast(overlaps);

        final JavaRDD<Shard<GATKRead>> readShards = createReadShards(shardBoundariesBroadcast, reads);

        final JavaRDD<Tuple2<AssemblyRegion, SimpleInterval>> assemblyRegions = readShards
                .mapPartitions(shardsToAssemblyRegions(authHolder, referenceBroadcast, hcArgsBroadcast, shardingArgs, header));

        return assemblyRegions.mapPartitions(callVariantsFromAssemblyRegions(authHolder, header, referenceBroadcast, hcArgsBroadcast));
    }

    /**
     * Call variants from Tuples of AssemblyRegion and Simple Interval
     * The interval should be the non-padded shard boundary for the shard that the corresponding AssemblyRegion was
     * created in, it's used to eliminate redundant variant calls at the edge of shard boundaries.
     */
    private static FlatMapFunction<Iterator<Tuple2<AssemblyRegion, SimpleInterval>>, VariantContext> callVariantsFromAssemblyRegions(
            final AuthHolder authHolder,
            final SAMFileHeader header,
            final Broadcast<ReferenceMultiSource> referenceBroadcast,
            final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast) {
        return regionAndIntervals -> {
            //HaplotypeCallerEngine isn't serializable but is expensive to instantiate, so construct and reuse one for every partition
            final ReferenceMultiSourceAdapter referenceReader = new ReferenceMultiSourceAdapter(referenceBroadcast.getValue(), authHolder);
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgsBroadcast.value(), header, referenceReader);
            return iteratorToStream(regionAndIntervals).flatMap(regionToVariants(hcEngine)).iterator();
        };
    }

    private static <T> Stream<T> iteratorToStream(Iterator<T> iterator) {
        Iterable<T> regionsIterable = () -> iterator;
        return StreamSupport.stream(regionsIterable.spliterator(), false);
    }

    private static Function<Tuple2<AssemblyRegion, SimpleInterval>, Stream<? extends VariantContext>> regionToVariants(HaplotypeCallerEngine hcEngine) {
        return regionAndInterval -> {
            final List<VariantContext> variantContexts = hcEngine.callRegion(regionAndInterval._1(), new FeatureContext());
            final SimpleInterval shardBoundary = regionAndInterval._2();
            return variantContexts.stream()
                .filter(vc -> shardBoundary.contains(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart())));
        };
    }

    /**
     * WriteVariants, this is currently going to be horribly slow and explosive on a full size file since it performs a collect.
     *
     * This will be replaced by a parallel writer similar to what's done with {@link org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink}
     */
    private void writeVariants(JavaRDD<VariantContext> variants) {
        final List<VariantContext> collectedVariants = variants.collect();
        final SAMSequenceDictionary referenceDictionary = getReferenceSequenceDictionary();

        final List<VariantContext> sortedVariants = collectedVariants.stream()
                .sorted((o1, o2) -> IntervalUtils.compareLocatables(o1, o2, referenceDictionary))
                .collect(Collectors.toList());

        final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, getHeaderForReads(), new ReferenceMultiSourceAdapter(getReference(), getAuthHolder()));
        try(final VariantContextWriter writer = hcEngine.makeVCFWriter(output, getBestAvailableSequenceDictionary())) {
            hcEngine.writeHeader(writer, getHeaderForReads().getSequenceDictionary(), Collections.emptySet());
            sortedVariants.forEach(writer::add);
        }
    }

    /**
     * Create an RDD of {@link Shard} from an RDD of {@link GATKRead}
     * @param shardBoundariesBroadcast  broadcast of an {@link OverlapDetector} loaded with the intervals that should be used for creating ReadShards
     * @param reads Rdd of {@link GATKRead}
     * @return a Rdd of reads grouped into potentially overlapping shards
     */
    private static JavaRDD<Shard<GATKRead>> createReadShards(final Broadcast<OverlapDetector<ShardBoundary>> shardBoundariesBroadcast, final JavaRDD<GATKRead> reads) {
        final JavaPairRDD<ShardBoundary, GATKRead> paired = reads.flatMapToPair(read -> {
            final Collection<ShardBoundary> overlappingShards = shardBoundariesBroadcast.value().getOverlaps(read);
            return overlappingShards.stream().map(key -> new Tuple2<>(key, read)).iterator();
        });
        final JavaPairRDD<ShardBoundary, Iterable<GATKRead>> shardsWithReads = paired.groupByKey();
        return shardsWithReads.map(shard -> new SparkReadShard(shard._1(), shard._2()));
    }

    /**
     * @return an {@link OverlapDetector} loaded with {@link ShardBoundary}
     * based on the -L intervals
     */
    private static OverlapDetector<ShardBoundary> getShardBoundaryOverlapDetector(final SAMFileHeader header, final List<SimpleInterval> intervals, final int readShardSize, final int readShardPadding) {
        final OverlapDetector<ShardBoundary> shardBoundaryOverlapDetector = new OverlapDetector<>(0, 0);
        intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, header.getSequenceDictionary()).stream())
                .forEach(boundary -> shardBoundaryOverlapDetector.addLhs(boundary, boundary.getPaddedInterval()));
        return shardBoundaryOverlapDetector;
    }

    /**
     * @return and RDD of {@link Tuple2<AssemblyRegion, SimpleInterval>} which pairs each AssemblyRegion with the
     * interval it was generated in
     */
    private static FlatMapFunction<Iterator<Shard<GATKRead>>, Tuple2<AssemblyRegion, SimpleInterval>> shardsToAssemblyRegions(
            final AuthHolder authHolder,
            final Broadcast<ReferenceMultiSource> reference,
            final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast,
            final ShardingArgumentCollection assemblyArgs,
            final SAMFileHeader header) {
        return shards -> {
            final ReferenceMultiSource referenceMultiSource = reference.value();
            final ReferenceMultiSourceAdapter referenceSource = new ReferenceMultiSourceAdapter(referenceMultiSource, authHolder);
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgsBroadcast.value(), header, referenceSource);

            return iteratorToStream(shards).flatMap(shardToRegion(assemblyArgs, header, referenceSource, hcEngine)).iterator();
        };
    }

    private static Function<Shard<GATKRead>, Stream<? extends Tuple2<AssemblyRegion, SimpleInterval>>> shardToRegion(
            ShardingArgumentCollection assemblyArgs,
            SAMFileHeader header,
            ReferenceMultiSourceAdapter referenceSource,
            HaplotypeCallerEngine evaluator) {
        return shard -> {
            final ReferenceContext refContext = new ReferenceContext(referenceSource, shard.getPaddedInterval());

            //TODO load features as a side input
            final FeatureContext features = new FeatureContext();

            final Iterable<AssemblyRegion> assemblyRegions = AssemblyRegion.createFromReadShard(
                    shard, header, refContext, features, evaluator,
                    assemblyArgs.minAssemblyRegionSize, assemblyArgs.maxAssemblyRegionSize,
                    assemblyArgs.assemblyRegionPadding, assemblyArgs.activeProbThreshold,
                    assemblyArgs.maxProbPropagationDistance);

            return StreamSupport.stream(assemblyRegions.spliterator(), false)
                    .map(a -> new Tuple2<>(a, shard.getInterval()));
        };
    }

    /**
     * Adapter to allow a 2bit reference to be used in HaplotypeCallerEngine.
     * This is not intended as a general purpose adapter, it only enables the operations needed in {@link HaplotypeCallerEngine}
     * This should not be used outside of this class except for testing purposes.
     */
    @VisibleForTesting
    public static final class ReferenceMultiSourceAdapter implements ReferenceSequenceFile, ReferenceDataSource, Serializable{
        private static final long serialVersionUID = 1L;

        private final ReferenceMultiSource source;
        private final AuthHolder auth;
        private final SAMSequenceDictionary sequenceDictionary;

        public ReferenceMultiSourceAdapter(final ReferenceMultiSource source, final AuthHolder auth) {
            this.source = source;
            this.auth = auth;
            sequenceDictionary = source.getReferenceSequenceDictionary(null);
        }

        @Override
        public ReferenceSequence queryAndPrefetch(final String contig, final long start, final long stop) {
           return getSubsequenceAt(contig, start, stop);
        }

        @Override
        public SAMSequenceDictionary getSequenceDictionary() {
            return source.getReferenceSequenceDictionary(null);
        }

        @Override
        public ReferenceSequence nextSequence() {
            throw new UnsupportedOperationException("nextSequence is not implemented");
        }

        @Override
        public void reset() {
            throw new UnsupportedOperationException("reset is not implemented");
        }

        @Override
        public boolean isIndexed() {
            return true;
        }

        @Override
        public ReferenceSequence getSequence(final String contig) {
            throw new UnsupportedOperationException("getSequence is not supported");
        }

        @Override
        public ReferenceSequence getSubsequenceAt(final String contig, final long start, final long stop) {
            try {
                final ReferenceBases bases = source.getReferenceBases(auth.asPipelineOptionsDeprecated(), new SimpleInterval(contig, (int) start, (int) stop));
                return new ReferenceSequence(contig, sequenceDictionary.getSequenceIndex(contig), bases.getBases());
            } catch (final IOException e) {
                throw new GATKException(String.format("Failed to load reference bases for %s:%d-%d", contig, start, stop));
            }
        }

        @Override
        public void close() {
            // doesn't do anything because you can't close a two-bit file
        }

        @Override
        public Iterator<Byte> iterator() {
            throw new UnsupportedOperationException("iterator is not supported");
        }
    }

}
