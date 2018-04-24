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
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.DuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import scala.Tuple2;

import javax.validation.constraints.Max;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@DocumentedFeature
@CommandLineProgramProperties(
        summary ="Marks duplicates on spark",
        oneLineSummary ="MarkDuplicates on Spark",
        programGroup = ReadDataManipulationProgramGroup.class)
@BetaFeature
public final class MarkDuplicatesSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    public static final String DO_NOT_MARK_UNMAPPED_MATES = "do-not-mark-unmapped-mates";

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

    @Argument(fullName = DO_NOT_MARK_UNMAPPED_MATES, doc = "Enabling this option will mean unmapped mates of duplicate marked reads will not be marked as duplicates.")
    public boolean dontMarkUnmappedMates = false;

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final MarkDuplicatesScoringStrategy scoringStrategy,
                                         final OpticalDuplicateFinder opticalDuplicateFinder,
                                         final int numReducers, final boolean dontMarkUnmappedMates) {

        JavaPairRDD<MarkDuplicatesSparkUtils.IndexPair<String>, Integer> namesOfNonDuplicates = MarkDuplicatesSparkUtils.transformToDuplicateNames(header, scoringStrategy, opticalDuplicateFinder, reads, numReducers);

        final JavaRDD<Tuple2<String,Integer>> repartitionedReadNames = namesOfNonDuplicates
                .mapToPair(pair -> new Tuple2<>(pair._1.getIndex(), new Tuple2<>(pair._1.getValue(),pair._2)))
                .partitionBy(new KnownIndexPartitioner(reads.getNumPartitions()))
                .values();
        
        return reads.zipPartitions(repartitionedReadNames, (readsIter, readNamesIter)  -> {
            final Map<String,Integer> namesOfNonDuplicateReadsAndOpticalCounts = Utils.stream(readNamesIter).collect(Collectors.toMap(Tuple2::_1,Tuple2::_2));
            return Utils.stream(readsIter).peek(read -> {
                // Handle reads that have been marked as non-duplicates (which also get tagged with optical duplicate summary statistics)
                if( namesOfNonDuplicateReadsAndOpticalCounts.containsKey(read.getName())) { //todo figure out if we should be marking the unmapped mates of duplicate reads as duplicates
                    read.setIsDuplicate(false);
                    if (!(dontMarkUnmappedMates && read.isUnmapped())) {
                        int dupCount = namesOfNonDuplicateReadsAndOpticalCounts.replace(read.getName(), -1);
                        if (dupCount > -1) {
                            ((SAMRecordToGATKReadAdapter) read).setTransientAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, dupCount);
                        }
                    }
                    // Mark unmapped read pairs as non-duplicates
                } else if (MarkDuplicatesSparkUtils.readAndMateAreUnmapped(read)) {
                    read.setIsDuplicate(false);
                    // Everything else is a duplicate
                } else{
                    if (!(dontMarkUnmappedMates && read.isUnmapped())) {
                        read.setIsDuplicate(true);
                    }
                }
            }).iterator();
        });
    }

    /**
     * A custom partitioner designed to cut down on spark shuffle costs.
     * This is designed such that getPartition(key) is called on a key which corresponds to the already known target partition
     */
    private static class KnownIndexPartitioner extends Partitioner {
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

        final SAMFileHeader header = getHeaderForReads();
        final JavaRDD<GATKRead> finalReadsForMetrics = mark(reads, header, duplicatesScoringStrategy, finder, getRecommendedNumReducers(), dontMarkUnmappedMates);

        if (metricsFile != null) {
            final JavaPairRDD<String, DuplicationMetrics> metricsByLibrary = MarkDuplicatesSparkUtils.generateMetrics(
                    header, finalReadsForMetrics);
            final MetricsFile<DuplicationMetrics, Double> resultMetrics = getMetricsFile();
            MarkDuplicatesSparkUtils.saveMetricsRDD(resultMetrics, header, metricsByLibrary, metricsFile);
        }
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        writeReads(ctx, output, finalReadsForMetrics, header);
    }

}
