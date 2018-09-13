package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import htsjdk.samtools.Defaults;
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
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.GATKDuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.SerializableOpticalDuplicatesFinder;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import scala.Tuple2;

import java.util.Collections;
import java.util.HashMap;
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

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "Path to write duplication metrics to.", optional=true,
            shortName = StandardArgumentDefinitions.METRICS_FILE_SHORT_NAME,
            fullName = StandardArgumentDefinitions.METRICS_FILE_LONG_NAME)
    protected String metricsFile;

    @ArgumentCollection
    protected MarkDuplicatesSparkArgumentCollection markDuplicatesSparkArgumentCollection = new MarkDuplicatesSparkArgumentCollection();

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    /**
     * Main method for marking duplicates, takes an JavaRDD of GATKRead and an associated SAMFileHeader with corresponding
     * sorting information and returns a new JavaRDD\<GATKRead\> in which all read templates have been marked as duplicates
     *
     * NOTE: This method expects the incoming reads to be grouped by read name (queryname sorted/querygrouped) and for this
     *       to be explicitly be set in the the provided header. Furthermore, all the reads in a template must be grouped
     *       into the same partition or there may be problems duplicate marking.
     *       If MarkDuplicates detects reads are sorted in some other way, it will perform an extra sort operation first,
     *       thus it is preferable to input reads to this method sorted for performance reasons.
     *
     * @param reads input reads to be duplicate marked
     * @param header header corresponding to the input reads
     * @param scoringStrategy method by which duplicates are detected
     * @param opticalDuplicateFinder
     * @param numReducers number of partitions to separate the data into
     * @param dontMarkUnmappedMates when true, unmapped mates of duplicate fragments will be marked as non-duplicates
     * @return A JavaRDD of GATKReads where duplicate flags have been set
     */
    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final MarkDuplicatesScoringStrategy scoringStrategy,
                                         final OpticalDuplicateFinder opticalDuplicateFinder,
                                         final int numReducers, final boolean dontMarkUnmappedMates) {
        final boolean markUnmappedMates = !dontMarkUnmappedMates;
        SAMFileHeader headerForTool = header.clone();

        // If the input isn't queryname sorted, sort it before duplicate marking
        final JavaRDD<GATKRead> sortedReadsForMarking = querynameSortReadsIfNecessary(reads, numReducers, headerForTool);

        final JavaPairRDD<MarkDuplicatesSparkUtils.IndexPair<String>, Integer> namesOfNonDuplicates = MarkDuplicatesSparkUtils.transformToDuplicateNames(headerForTool, scoringStrategy, opticalDuplicateFinder, sortedReadsForMarking, numReducers);

        // Here we explicitly repartition the read names of the unmarked reads to match the partitioning of the original bam
        final JavaRDD<Tuple2<String,Integer>> repartitionedReadNames = namesOfNonDuplicates
                .mapToPair(pair -> new Tuple2<>(pair._1.getIndex(), new Tuple2<>(pair._1.getValue(),pair._2)))
                .partitionBy(new KnownIndexPartitioner(sortedReadsForMarking.getNumPartitions()))
                .values();

        // Here we combine the original bam with the repartitioned unmarked readnames to produce our marked reads
        return sortedReadsForMarking.zipPartitions(repartitionedReadNames, (readsIter, readNamesIter)  -> {
            final Map<String,Integer> namesOfNonDuplicateReadsAndOpticalCounts = new HashMap<>();
            readNamesIter.forEachRemaining(tup -> { if (namesOfNonDuplicateReadsAndOpticalCounts.putIfAbsent(tup._1,tup._2)!=null) {
                throw new GATKException(String.format("Detected multiple mark duplicate records objects corresponding to read with name '%s', this could be the result of the file sort order being incorrect or that a previous tool has let readnames span multiple partitions",tup._1()));
            }
                });

            return Utils.stream(readsIter)
                    .peek(read -> read.setIsDuplicate(false))
                    .peek(read -> {
                // Handle reads that have been marked as non-duplicates (which also get tagged with optical duplicate summary statistics)
                if (namesOfNonDuplicateReadsAndOpticalCounts.containsKey(read.getName())) {
                    read.setIsDuplicate(false);
                    if (markUnmappedMates || !read.isUnmapped()) {
                        int dupCount = namesOfNonDuplicateReadsAndOpticalCounts.replace(read.getName(), -1);
                        if (dupCount > -1) {
                            ((SAMRecordToGATKReadAdapter) read).setTransientAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, dupCount);
                        }
                    }
                    // Mark unmapped read pairs as non-duplicates
                } else if (ReadUtils.readAndMateAreUnmapped(read)) {
                    read.setIsDuplicate(false);
                    // Everything else is a duplicate
                } else {
                    if (markUnmappedMates || !read.isUnmapped()) {
                        read.setIsDuplicate(true);
                    } else {
                        read.setIsDuplicate(false);
                    }
                }
            }).iterator();
        });
    }

    /**
     * Sort reads into queryname order if they are not already sorted
     */
    private static JavaRDD<GATKRead> querynameSortReadsIfNecessary(JavaRDD<GATKRead> reads, int numReducers, SAMFileHeader headerForTool) {
        JavaRDD<GATKRead> sortedReadsForMarking;
        if (ReadUtils.isReadNameGroupedBam(headerForTool)) {
            sortedReadsForMarking = reads;
        } else {
            headerForTool.setSortOrder(SAMFileHeader.SortOrder.queryname);
            sortedReadsForMarking = SparkUtils.sortReadsAccordingToHeader(reads, headerForTool, numReducers);
        }
        return sortedReadsForMarking;
    }

    /**
     * A custom partitioner designed to cut down on spark shuffle costs.
     * This is designed such that getPartition(key) is called on a key which corresponds to the already known target partition
     *
     * By storing the original partitioning for each read and passing it through the duplicates marking process it
     * allows us to get away with just shuffling the small read name objects to the correct partition in the original bam
     * while avoiding any shuffle of the larger read objects.
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
                new SerializableOpticalDuplicatesFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE) : null;

        final SAMFileHeader header = getHeaderForReads();
        final JavaRDD<GATKRead> finalReadsForMetrics = mark(reads, header, markDuplicatesSparkArgumentCollection.duplicatesScoringStrategy, finder, getRecommendedNumReducers(),  markDuplicatesSparkArgumentCollection.dontMarkUnmappedMates);

        if (metricsFile != null) {
            final JavaPairRDD<String, GATKDuplicationMetrics> metricsByLibrary = MarkDuplicatesSparkUtils.generateMetrics(
                    header, finalReadsForMetrics);
            final MetricsFile<GATKDuplicationMetrics, Double> resultMetrics = getMetricsFile();
            MarkDuplicatesSparkUtils.saveMetricsRDD(resultMetrics, header, metricsByLibrary, metricsFile);
        }
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        writeReads(ctx, output, finalReadsForMetrics, header);
    }

}
