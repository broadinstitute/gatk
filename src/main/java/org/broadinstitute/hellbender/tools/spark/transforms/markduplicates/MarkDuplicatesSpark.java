package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import scala.Tuple2;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@DocumentedFeature
@CommandLineProgramProperties(
        summary ="Marks duplicates on spark",
        oneLineSummary ="MarkDuplicates on Spark",
        programGroup = ReadDataManipulationProgramGroup.class)
@BetaFeature
public final class MarkDuplicatesSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "Path to write duplication metrics to.", optional=true,
            shortName = "M", fullName = "METRICS_FILE")
    protected String metricsFile;

    @Argument(shortName = "DS", fullName = "DUPLICATE_SCORING_STRATEGY", doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final MarkDuplicatesScoringStrategy scoringStrategy,
                                         final OpticalDuplicateFinder opticalDuplicateFinder, final int numReducers) {

     //   JavaRDD<GATKRead> primaryReads = reads.filter(v1 -> !ReadUtils.isNonPrimary(v1));
//        JavaRDD<GATKRead> nonPrimaryReads = reads.filter(v1 -> ReadUtils.isNonPrimary(v1));
        JavaPairRDD<MarkDuplicatesSparkUtils.IndexPair<String>, Integer> namesOfNonDuplicates = MarkDuplicatesSparkUtils.transformToDuplicateNames(header, scoringStrategy, opticalDuplicateFinder, reads, numReducers);

        final JavaRDD<String> repartitionedReadNames = namesOfNonDuplicates.keys()
                .mapToPair(pair -> new Tuple2<>(pair.getIndex(), pair.getValue()))
                .partitionBy(new KnownIndexPartitioner(reads.getNumPartitions()))
                .values();
        
        return reads.zipPartitions(repartitionedReadNames, (readsIter, readNamesIter)  -> {
            final Set<String> namesOfNonDuplicateReads = Utils.stream(readNamesIter).collect(Collectors.toSet());
            return Utils.stream(readsIter).map( read -> {
                if( namesOfNonDuplicateReads.contains(read.getName())) {
                    read.setIsDuplicate(false);
                } else {
                    read.setIsDuplicate(true);
                }
                return read;
            }).iterator();
        });
    }

    public static class KnownIndexPartitioner extends Partitioner {
        private static final long serialVersionUID = 1L;
        private final int numPartitions;

        KnownIndexPartitioner(int numPartitions) {
            this.numPartitions = numPartitions;
        }

        @Override
        public int numPartitions() {
            return numPartitions;
        }

        @Override
        @SuppressWarnings("unchecked")
        public int getPartition(Object key) {
            return (Integer) key;
        }
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> reads = getReads();
        final OpticalDuplicateFinder finder = opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null ?
                new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;
        //reads.partitions().get(0).
        //.mapToPair(read -> new Tuple2<String, GATKRead>(ReadsKey.keyForRead(getHeaderForReads(), read), read)

        final JavaRDD<GATKRead> finalReadsForMetrics = mark(reads, getHeaderForReads(), duplicatesScoringStrategy, finder, getRecommendedNumReducers());
        //finalReadsForMetrics.repartition(new RangePartitioner<>())
        if (metricsFile != null) {
            final JavaPairRDD<String, DuplicationMetrics> metricsByLibrary = MarkDuplicatesSparkUtils.generateMetrics(getHeaderForReads(), finalReadsForMetrics);
            final MetricsFile<DuplicationMetrics, Double> resultMetrics = getMetricsFile();
            MarkDuplicatesSparkUtils.saveMetricsRDD(resultMetrics, getHeaderForReads(), metricsByLibrary, metricsFile);
        }

        final JavaRDD<GATKRead> finalReads = cleanupTemporaryAttributes( finalReadsForMetrics);
        writeReads(ctx, output, finalReads);
    }


    /**
     * The OD attribute was added to each read for optical dups.
     * Now we have to clear it to avoid polluting the output.
     */
    public static JavaRDD<GATKRead> cleanupTemporaryAttributes(final JavaRDD<GATKRead> reads) {
        return reads.map(read -> {
            read.clearAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME);
            return read;
        });
    }
}
