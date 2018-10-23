package org.broadinstitute.hellbender.tools;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.spark.SparkFiles;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionReadShardArgumentCollection;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.ShardToMultiIntervalShardAdapter;
import org.broadinstitute.hellbender.engine.spark.SparkSharder;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.PositionalDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import scala.Tuple2;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * ********************************************************************************
 * *   This tool DOES NOT match the output of HaplotypeCaller.                    *
 * *   It is still under development and should not be used for production work.  *
 * *   For evaluation only.                                                       *
 * *   Use the non-spark HaplotypeCaller if you care about the results.           *
 * ********************************************************************************
 *
 * Call germline SNPs and indels via local re-assembly of haplotypes.
 *
 * <p>This is an implementation of {@link HaplotypeCaller} using spark to distribute the computation.
 * It is still in an early stage of development and does not yet support all the options that the non-spark version does.
 * Specifically it does not support the --dbsnp, --comp, and --bam-output options.</p>
 *
 * <h3>Usage Example</h3>
 * <pre>
 * gatk HaplotypeCallerSpark \
 * -R Homo_sapiens_assembly38.fasta \
 * -I input.bam \
 * -O output.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(summary = "HaplotypeCaller on Spark", oneLineSummary = "HaplotypeCaller on Spark", programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
@BetaFeature
public final class HaplotypeCallerSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    public static final int DEFAULT_READSHARD_SIZE = 5000;
    private static final boolean INCLUDE_READS_WITH_DELETIONS_IN_IS_ACTIVE_PILEUPS = true;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Single file to which variants should be written")
    public String output;

    public static class HaplotypeCallerAssemblyRegionArgumentCollection extends AssemblyRegionArgumentCollection {
        private static final long serialVersionUID = 1L;

        @Override
        protected int defaultMinAssemblyRegionSize() { return HaplotypeCaller.DEFAULT_MIN_ASSEMBLY_REGION_SIZE; }

        @Override
        protected int defaultMaxAssemblyRegionSize() { return HaplotypeCaller.DEFAULT_MAX_ASSEMBLY_REGION_SIZE; }

        @Override
        protected int defaultAssemblyRegionPadding() { return HaplotypeCaller.DEFAULT_ASSEMBLY_REGION_PADDING; }

        @Override
        protected int defaultMaxReadsPerAlignmentStart() { return HaplotypeCaller.DEFAULT_MAX_READS_PER_ALIGNMENT; }

        @Override
        protected double defaultActiveProbThreshold() { return HaplotypeCaller.DEFAULT_ACTIVE_PROB_THRESHOLD; }

        @Override
        protected int defaultMaxProbPropagationDistance() { return HaplotypeCaller.DEFAULT_MAX_PROB_PROPAGATION_DISTANCE; }
    }

    @ArgumentCollection
    public final AssemblyRegionReadShardArgumentCollection shardingArgs = new AssemblyRegionReadShardArgumentCollection();

    @ArgumentCollection
    public final AssemblyRegionArgumentCollection assemblyRegionArgs = new HaplotypeCallerAssemblyRegionArgumentCollection();

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
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() { return HaplotypeCallerEngine.getStandardHaplotypeCallerAnnotationGroups();}

    @Override
    public Collection<Annotation> makeVariantAnnotations() {
        final boolean referenceConfidenceMode = hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE;
        final Collection<Annotation> annotations = super.makeVariantAnnotations();
        return referenceConfidenceMode? HaplotypeCallerEngine.filterReferenceConfidenceAnnotations(annotations): annotations;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        //TODO remove me when https://github.com/broadinstitute/gatk/issues/4303 are fixed
        if (output.endsWith(IOUtil.BCF_FILE_EXTENSION) || output.endsWith(IOUtil.BCF_FILE_EXTENSION + ".gz")) {
            throw new UserException.UnimplementedFeature("It is currently not possible to write a BCF file on spark.  See https://github.com/broadinstitute/gatk/issues/4303 for more details .");
        }
        SAMFileHeader header = getHeaderForReads();
        Utils.validate(header.getSortOrder() == SAMFileHeader.SortOrder.coordinate, "The reads must be coordinate sorted.");
        logger.info("********************************************************************************");
        logger.info("The output of this tool DOES NOT match the output of HaplotypeCaller. ");
        logger.info("It is under development and should not be used for production work. ");
        logger.info("For evaluation only.");
        logger.info("Use the non-spark HaplotypeCaller if you care about the results. ");
        logger.info("********************************************************************************");
        addReferenceFilesForSpark(ctx, referenceArguments.getReferenceFileName());
        final List<SimpleInterval> intervals = hasUserSuppliedIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(header.getSequenceDictionary());
        callVariantsWithHaplotypeCallerAndWriteOutput(ctx, getReads(), header, referenceArguments.getReferenceFileName(), intervals, hcArgs, shardingArgs, assemblyRegionArgs, numReducers, output, makeVariantAnnotations());
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return HaplotypeCallerEngine.makeStandardHCReadFilters();
    }

    /**
     * Call Variants using HaplotypeCaller on Spark and write out a VCF file.
     *
     * This may be called from any spark pipeline in order to call variants from an RDD of GATKRead
     *
     * @param ctx the spark context
     * @param reads the reads variants should be called from
     * @param header the header that goes with the reads
     * @param reference the full path to the reference file (must have been added via {@code SparkContext#addFile()})
     * @param intervals the intervals to restrict calling to
     * @param hcArgs haplotype caller arguments
     * @param shardingArgs arguments to control how the assembly regions are sharded
     * @param numReducers the number of reducers to use when sorting
     * @param output the output path for the VCF
     */
    public static void callVariantsWithHaplotypeCallerAndWriteOutput(
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads,
            final SAMFileHeader header,
            final String reference,
            final List<SimpleInterval> intervals,
            final HaplotypeCallerArgumentCollection hcArgs,
            final AssemblyRegionReadShardArgumentCollection shardingArgs,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final int numReducers,
            final String output,
            final Collection<Annotation> annotations) {
        final VariantAnnotatorEngine variantannotatorEngine = new VariantAnnotatorEngine(annotations,  hcArgs.dbsnp.dbsnp, hcArgs.comps, hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE);

        final Path referencePath = IOUtils.getPath(reference);
        final ReferenceSequenceFile driverReferenceSequenceFile = new CachingIndexedFastaSequenceFile(referencePath);
        final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, false, false, header, driverReferenceSequenceFile, variantannotatorEngine);
        final String referenceFileName = referencePath.getFileName().toString();
        final JavaRDD<VariantContext> variants = callVariantsWithHaplotypeCaller(ctx, reads, header, referenceFileName, intervals, hcArgs, shardingArgs, assemblyRegionArgs, variantannotatorEngine);
        variants.cache(); // without caching, computations are run twice as a side effect of finding partition boundaries for sorting
        try {
            VariantsSparkSink.writeVariants(ctx, output, variants, hcEngine.makeVCFHeader(header.getSequenceDictionary(), new HashSet<>()),
                    hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF, hcArgs.GVCFGQBands, hcArgs.genotypeArgs.samplePloidy);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, "writing failed", e);
        }
    }

    /**
     * Call Variants using HaplotypeCaller on Spark and return an RDD of  {@link VariantContext}
     *
     * This may be called from any spark pipeline in order to call variants from an RDD of GATKRead
     *
     * @param ctx the spark context
     * @param reads the reads variants should be called from
     * @param header the header that goes with the reads
     * @param referenceFileName the name of the reference file added via {@code SparkContext#addFile()}
     * @param intervals the intervals to restrict calling to
     * @param hcArgs haplotype caller arguments
     * @param shardingArgs arguments to control how the assembly regions are sharded
     * @param variantannotatorEngine
     * @return an RDD of Variants
     */
    public static JavaRDD<VariantContext> callVariantsWithHaplotypeCaller(
            final JavaSparkContext ctx,
            final JavaRDD<GATKRead> reads,
            final SAMFileHeader header,
            final String referenceFileName,
            final List<SimpleInterval> intervals,
            final HaplotypeCallerArgumentCollection hcArgs,
            final AssemblyRegionReadShardArgumentCollection shardingArgs,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final VariantAnnotatorEngine variantannotatorEngine) {
        Utils.validateArg(hcArgs.dbsnp.dbsnp == null, "HaplotypeCallerSpark does not yet support -D or --dbsnp arguments" );
        Utils.validateArg(hcArgs.comps.isEmpty(), "HaplotypeCallerSpark does not yet support -comp or --comp arguments" );
        Utils.validateArg(hcArgs.bamOutputPath == null, "HaplotypeCallerSpark does not yet support -bamout or --bamOutput");

        final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast = ctx.broadcast(hcArgs);

        final Broadcast<VariantAnnotatorEngine> annotatorEngineBroadcast = ctx.broadcast(variantannotatorEngine);

        final List<ShardBoundary> shardBoundaries = getShardBoundaries(header, intervals, shardingArgs.readShardSize, shardingArgs.readShardPadding);

        final int maxReadLength = reads.map(r -> r.getEnd() - r.getStart() + 1).reduce(Math::max);

        final JavaRDD<Shard<GATKRead>> readShards = SparkSharder.shard(ctx, reads, GATKRead.class, header.getSequenceDictionary(), shardBoundaries, maxReadLength);

        final JavaRDD<Tuple2<AssemblyRegion, SimpleInterval>> assemblyRegions = readShards
                .mapPartitions(shardsToAssemblyRegions(referenceFileName,
                                                       hcArgsBroadcast, assemblyRegionArgs, header, annotatorEngineBroadcast));

        return assemblyRegions.mapPartitions(callVariantsFromAssemblyRegions(header, referenceFileName, hcArgsBroadcast, annotatorEngineBroadcast));
    }

    /**
     * Call variants from Tuples of AssemblyRegion and Simple Interval
     * The interval should be the non-padded shard boundary for the shard that the corresponding AssemblyRegion was
     * created in, it's used to eliminate redundant variant calls at the edge of shard boundaries.
     */
    private static FlatMapFunction<Iterator<Tuple2<AssemblyRegion, SimpleInterval>>, VariantContext> callVariantsFromAssemblyRegions(
            final SAMFileHeader header,
            final String referenceFileName,
            final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast,
            final Broadcast<VariantAnnotatorEngine> annotatorEngineBroadcast) {
        return regionAndIntervals -> {
            //HaplotypeCallerEngine isn't serializable but is expensive to instantiate, so construct and reuse one for every partition
            final String pathOnExecutor = SparkFiles.get(referenceFileName);
            final ReferenceSequenceFile taskReferenceSequenceFile = new CachingIndexedFastaSequenceFile(IOUtils.getPath(pathOnExecutor));
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgsBroadcast.value(), false, false, header, taskReferenceSequenceFile, annotatorEngineBroadcast.getValue());
            return Utils.stream(regionAndIntervals).flatMap(regionToVariants(hcEngine)).iterator();
        };
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
     * @return a list of {@link ShardBoundary}
     * based on the -L intervals
     */
    private static List<ShardBoundary> getShardBoundaries(final SAMFileHeader
        header, final List<SimpleInterval> intervals, final int readShardSize, final int readShardPadding) {
        return intervals.stream()
            .flatMap(interval -> Shard.divideIntervalIntoShards(interval, readShardSize, readShardPadding, header.getSequenceDictionary()).stream())
            .collect(Collectors.toList());
    }

    /**
     * @return and RDD of {@link Tuple2<AssemblyRegion, SimpleInterval>} which pairs each AssemblyRegion with the
     * interval it was generated in
     */
    private static FlatMapFunction<Iterator<Shard<GATKRead>>, Tuple2<AssemblyRegion, SimpleInterval>> shardsToAssemblyRegions(
            final String referenceFileName,
            final Broadcast<HaplotypeCallerArgumentCollection> hcArgsBroadcast,
            final AssemblyRegionArgumentCollection assemblyRegionArgs,
            final SAMFileHeader header,
            final Broadcast<VariantAnnotatorEngine> annotatorEngineBroadcast) {
        return shards -> {
            final String pathOnExecutor = SparkFiles.get(referenceFileName);
            final ReferenceSequenceFile taskReferenceSequenceFile = new CachingIndexedFastaSequenceFile(IOUtils.getPath(pathOnExecutor));
            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgsBroadcast.value(), false, false, header, taskReferenceSequenceFile, annotatorEngineBroadcast.getValue());

            final ReferenceDataSource taskReferenceDataSource = new ReferenceFileSource(IOUtils.getPath(pathOnExecutor)); // TODO: share with taskReferenceSequenceFile
            final ReadsDownsampler readsDownsampler = assemblyRegionArgs.maxReadsPerAlignmentStart > 0 ?
                new PositionalDownsampler(assemblyRegionArgs.maxReadsPerAlignmentStart, header) : null;
            return Utils.stream(shards)
                    //TODO we've hacked multi interval shards here with a shim, but we should investigate as smarter approach https://github.com/broadinstitute/gatk/issues/4299
                .map(shard -> new ShardToMultiIntervalShardAdapter<>(
                        new DownsampleableSparkReadShard(new ShardBoundary(shard.getInterval(), shard.getPaddedInterval()), shard, readsDownsampler)))
                .flatMap(shardToRegion(assemblyRegionArgs, header, taskReferenceDataSource, hcEngine)).iterator();
        };
    }

    private static Function<MultiIntervalShard<GATKRead>, Stream<? extends Tuple2<AssemblyRegion, SimpleInterval>>> shardToRegion(
            AssemblyRegionArgumentCollection assemblyRegionArgs,
            SAMFileHeader header,
            ReferenceDataSource referenceDataSource,
            HaplotypeCallerEngine evaluator) {
        return shard -> {
            //TODO load features as a side input
            final FeatureManager featureManager = null;

            final Iterator<AssemblyRegion> assemblyRegionIter = new AssemblyRegionIterator(shard, header, referenceDataSource, featureManager, evaluator, assemblyRegionArgs.minAssemblyRegionSize, assemblyRegionArgs.maxAssemblyRegionSize, assemblyRegionArgs.assemblyRegionPadding, assemblyRegionArgs.activeProbThreshold, assemblyRegionArgs.maxProbPropagationDistance,
                                                                                           INCLUDE_READS_WITH_DELETIONS_IN_IS_ACTIVE_PILEUPS);

            return Utils.stream(assemblyRegionIter)
                    .map(a -> new Tuple2<>(a, shard.getIntervals().get(0)));
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

        private final ReferenceMultiSparkSource source;
        private final SAMSequenceDictionary sequenceDictionary;

        public ReferenceMultiSourceAdapter(final ReferenceMultiSparkSource source) {
            this.source = source;
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
                final ReferenceBases bases = source.getReferenceBases(new SimpleInterval(contig, (int) start, (int) stop));
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
