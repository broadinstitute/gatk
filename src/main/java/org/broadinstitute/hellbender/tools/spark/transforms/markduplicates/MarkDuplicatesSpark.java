package org.broadinstitute.hellbender.tools.spark.transforms.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.Partitioner;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.GATKDuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.markduplicates.MarkDuplicates;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import scala.Tuple2;

import java.util.*;

/**
 * MarkDuplicates on Spark
 *
 * <p>This is a Spark implementation of <a href='https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php'>Picard MarkDuplicates</a> that allows the tool to be run in parallel on multiple cores on a local machine or multiple machines on a Spark cluster while still matching the output of the non-Spark Picard version of the tool. Since the tool requires holding all of the readnames in memory while it groups read information, machine configuration and starting sort-order impact tool performance. </p>
 *
 * Here are some differences of note between MarkDuplicatesSpark and Picard MarkDuplicates.
 *
 * <ul>
 *  <li>MarkDuplicatesSpark processing can replace both the MarkDuplicates and SortSam steps of the Best Practices <a href="https://software.broadinstitute.org/gatk/documentation/article?id=7899#2">single sample pipeline</a>. After flagging duplicate sets, the tool automatically coordinate-sorts the records. It is still necessary to subsequently run SetNmMdAndUqTags before running BQSR. </li>
 *  <li>The tool is optimized to run on queryname-grouped alignments (that is, all reads with the same queryname are together in the input file). If provided coordinate-sorted alignments, the tool will spend additional time first queryname sorting the reads internally. This can result in the tool being up to 2x slower processing under some circumstances.</li>
 *  <li>Due to MarkDuplicatesSpark queryname-sorting coordinate-sorted inputs internally at the start, the tool produces identical results regardless of the input sort-order. That is, it will flag duplicates sets that include secondary, and supplementary and unmapped mate records no matter the sort-order of the input. This differs from how Picard MarkDuplicates behaves given the differently sorted inputs. </li>
 *  <li>Collecting duplicate metrics slows down performance and thus the metrics collection is optional and must be specified for the Spark version of the tool with '-M'. It is possible to collect the metrics with the standalone Picard tool <a href='https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_EstimateLibraryComplexity.php'>EstimateLibraryComplexity</a>.</li>
 *  <li>MarkDuplicatesSpark is optimized to run locally on a single machine by leveraging core parallelism that MarkDuplicates and SortSam cannot. It will typically run faster than MarkDuplicates and SortSam by a factor of 15% over the same data at 2 cores and will scale linearly to upwards of 16 cores. This means MarkDuplicatesSpark, even without access to a Spark cluster, is faster than MarkDuplicates.</li>
 *  <li>MarkDuplicatesSpark can be run with multiple input bams. If this is the case all of the inputs must be a mix queryname-grouped or queryname sorted.</li>
 * </ul>
 *
 * <p>For a typical 30x coverage WGS BAM, we recommend running on a machine with at least 16 GB. Memory usage scales with library complexity and the tool will need more memory for larger or more complex data. If the tool is running slowly it is possible Spark is running out of memory and is spilling data to disk excessively. If this is the case then increasing the memory available to the tool should yield speedup to a threshold; otherwise, increasing memory should have no effect beyond that threshold. </p>
 *
 * <p> Note that this tool does not support UMI based duplicate marking. </p>
 *
 * <p>See <a href='https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php'>MarkDuplicates documentation</a> for details on tool features and background information. </p>
 *
 * <h3>Usage examples</h3>
 * Provide queryname-grouped reads to MarkDuplicatesSpark
 *     <pre>
 *      gatk MarkDuplicatesSpark \
 *            -I input.bam \
 *            -O marked_duplicates.bam
 *     </pre>
 *
 * Additionally produce estimated library complexity metrics
 *     <pre>
 *     gatk MarkDuplicatesSpark \
 *             -I input.bam \
 *             -O marked_duplicates.bam \
 *             -M marked_dup_metrics.txt
 *
 *     </pre>
 *
 *
 * MarkDuplicatesSpark run locally specifying the removal of sequencing duplicates
 *     <pre>
 *       gatk MarkDuplicatesSpark \
 *            -I input.bam \
 *            -O marked_duplicates.bam \
 *            --remove-sequencing-duplicates
 *     </pre>
 *
 * MarkDuplicatesSpark run locally tagging OpticalDuplicates using the "DT" attribute for reads
 *     <pre>
 *       gatk MarkDuplicatesSpark \
 *            -I input.bam \
 *            -O marked_duplicates.bam \
 *            --duplicate-tagging-policy OpticalOnly
 *     </pre>
 *
 *  MarkDuplicates run locally specifying the core input. Note if 'spark.executor.cores' is unset, Spark will use all available cores on the machine.
 *     <pre>
 *       gatk MarkDuplicatesSpark \
 *            -I input.bam \
 *            -O marked_duplicates.bam \
 *            -M marked_dup_metrics.txt \
 *            --conf 'spark.executor.cores=5'
 *     </pre>
 *
 *  MarkDuplicates run on a Spark cluster of five executors  and with eight executor cores
 *     <pre>
 *       gatk MarkDuplicatesSpark \
 *            -I input.bam \
 *            -O marked_duplicates.bam \
 *            -M marked_dup_metrics.txt \
 *            -- \
 *            --spark-runner SPARK \
 *            --spark-master MASTER_URL \
 *            --num-executors 5 \
 *            --executor-cores 8
 *     </pre>
 *
 *    Please see
 *    <a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics'>Picard DuplicationMetrics</a>
 *    for detailed explanations of the output metrics.
 *    <hr />
 *
 * <h3>Notes</h3>
 * <ol>
 *     <li>This Spark tool requires a significant amount of disk operations. Run with both the input data and outputs on high throughput SSDs when possible. When pipelining this tool on Google Compute Engine instances, for best performance requisition machines with LOCAL SSDs.  </li>
 *     <li>Furthermore, we recommend explicitly setting the Spark temp directory to an available SSD when running this in local mode by adding the argument --conf 'spark.local.dir=/PATH/TO/TEMP/DIR'. See <a href="https://gatkforums.broadinstitute.org/gatk/discussion/comment/56337">this forum discussion</a> for details.</li>
 * </ol>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary ="Marks duplicates on spark",
        oneLineSummary ="MarkDuplicates on Spark",
        programGroup = ReadDataManipulationProgramGroup.class)
public final class MarkDuplicatesSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    public static final String ALLOW_MULTIPLE_SORT_ORDERS_IN_INPUT_ARG = "allow-multiple-sort-orders-in-input";
    public static final String TREAT_UNSORTED_AS_ORDERED = "treat-unsorted-as-querygroup-ordered";

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "Path to write duplication metrics to.", optional=true,
            shortName = StandardArgumentDefinitions.METRICS_FILE_SHORT_NAME,
            fullName = StandardArgumentDefinitions.METRICS_FILE_LONG_NAME)
    protected String metricsFile;

    @Advanced
    @Argument(doc = "Allow non-queryname sorted inputs when specifying multiple input bams.", optional=true,
            fullName = ALLOW_MULTIPLE_SORT_ORDERS_IN_INPUT_ARG)
    protected boolean allowMultipleSortOrders = false;

    @Advanced
    @Argument(doc = "Treat unsorted files as query-group orderd files. WARNING: This option disables a basic safety check and may result in unexpected behavior if the file is truly unordered", optional=true,
            fullName = TREAT_UNSORTED_AS_ORDERED)
    protected boolean treatUnsortedAsOrdered = false;

    @ArgumentCollection
    protected MarkDuplicatesSparkArgumentCollection markDuplicatesSparkArgumentCollection = new MarkDuplicatesSparkArgumentCollection();

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    // Reads with this marker will be treated as non-duplicates always
    public static int NO_OPTICAL_MARKER = -1;
    // Reads with this marker will be treated and marked as optical duplicates
    public static int OPTICAL_DUPLICATE_MARKER = -2;

    @Override
    public ReadInputMergingPolicy getReadInputMergingPolicy() {
        return ReadInputMergingPolicy.concatMerge;
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
     * @param taggingPolicy determines whether optical duplicates and library duplicates are labeled with the "DT" tag
     * @return A JavaRDD of GATKReads where duplicate flags have been set
     */
    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final MarkDuplicatesScoringStrategy scoringStrategy,
                                         final OpticalDuplicateFinder opticalDuplicateFinder,
                                         final int numReducers, final boolean dontMarkUnmappedMates,
                                         final MarkDuplicates.DuplicateTaggingPolicy taggingPolicy) {
        final boolean markUnmappedMates = !dontMarkUnmappedMates;
        SAMFileHeader headerForTool = header.clone();

        // If the input isn't queryname sorted, sort it before duplicate marking
        final JavaRDD<GATKRead> sortedReadsForMarking = SparkUtils.querynameSortReadsIfNecessary(reads, numReducers, headerForTool);

        // If we need to remove optical duplicates or tag them, then make sure we are keeping track
        final boolean markOpticalDups = (taggingPolicy != MarkDuplicates.DuplicateTaggingPolicy.DontTag);

        final JavaPairRDD<MarkDuplicatesSparkUtils.IndexPair<String>, Integer> namesOfNonDuplicates = MarkDuplicatesSparkUtils.transformToDuplicateNames(headerForTool, scoringStrategy, opticalDuplicateFinder, sortedReadsForMarking, numReducers, markOpticalDups);

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
                    .peek(read -> read.setAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG, (String) null))
                    .peek(read -> {
                        // Handle reads that have been marked as non-duplicates (which also get tagged with optical duplicate summary statistics)
                        if (namesOfNonDuplicateReadsAndOpticalCounts.containsKey(read.getName())) {
                            // If its an optical duplicate, mark it. (Note: we only expect these to exist if optical duplicate marking is on)
                            if (namesOfNonDuplicateReadsAndOpticalCounts.get(read.getName()) == OPTICAL_DUPLICATE_MARKER) {
                                read.setIsDuplicate(true);
                                read.setAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG, MarkDuplicates.DUPLICATE_TYPE_SEQUENCING);

                            // Otherwise treat it normally as a non-duplicate.
                            } else {
                                read.setIsDuplicate(false);
                                if (markUnmappedMates || !read.isUnmapped()) {
                                    int dupCount = namesOfNonDuplicateReadsAndOpticalCounts.replace(read.getName(), NO_OPTICAL_MARKER);
                                    if (dupCount > -1) {
                                        read.setTransientAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, dupCount);
                                    }
                                }
                            }
                            // Mark unmapped read pairs as non-duplicates
                        } else if (ReadUtils.readAndMateAreUnmapped(read)) {
                            read.setIsDuplicate(false);
                            // Everything else is a duplicate
                        } else {
                            if (markUnmappedMates || !read.isUnmapped()) {
                                read.setIsDuplicate(true);
                                if (taggingPolicy == MarkDuplicates.DuplicateTaggingPolicy.All) {
                                    read.setAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG, MarkDuplicates.DUPLICATE_TYPE_LIBRARY);
                                }
                            } else {
                                read.setIsDuplicate(false);
                            }
                        }
                    }).iterator();
        });
    }

    public static JavaRDD<GATKRead> mark(final JavaRDD<GATKRead> reads, final SAMFileHeader header,
                                         final OpticalDuplicateFinder finder,
                                         final MarkDuplicatesSparkArgumentCollection mdArgs,
                                         final int numReducers) {
        return mark(reads,
                    header,
                    mdArgs.duplicatesScoringStrategy,
                    finder,
                    numReducers,
                    mdArgs.dontMarkUnmappedMates,
                    mdArgs.taggingPolicy);
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
        final SAMFileHeader mergedHeader = getHeaderForReads();

        // Check if we are using multiple inputs that the headers are all in the correct querygrouped ordering, if so set the aggregate header to reflect this
        if (readArguments.getReadPathSpecifiers().size() > 1) {
            final Optional<GATKPath> badlySorted = readArguments.getReadPathSpecifiers().stream()
                    .filter(spec -> !treatAsReadGroupOrdered(getHeaderForReadsInput(spec), treatUnsortedAsOrdered))
                    .findFirst();

            if(badlySorted.isPresent()) {
                if (allowMultipleSortOrders) {
                    //don't set an ordering, the files will all be sorted downstream
                    logger.info("Input files are not all grouped by read name so they will be sorted together.");
                } else {
                    throw new UserException(
                            "Multiple inputs to MarkDuplicatesSpark detected. MarkDuplicatesSpark requires all inputs to be queryname sorted " +
                                    "or querygroup-sorted for multi-input processing but input " + badlySorted.get() + " was sorted in " +
                                    getHeaderForReadsInput(badlySorted.get()) + " order");
                }
            } else {
                // The default sort order for merged input files is unsorted, so this will be fed to the tool to be sorted
                if (!allowMultipleSortOrders) {
                    mergedHeader.setGroupOrder(SAMFileHeader.GroupOrder.query);
                }
            }

        // If there is only one file and we are in treatUnsortedAsOrdered mode than set its group order accordingly.
        } else {
            if (treatUnsortedAsOrdered && (mergedHeader.getSortOrder().equals(SAMFileHeader.SortOrder.unknown) || mergedHeader.getSortOrder().equals(SAMFileHeader.SortOrder.unsorted))) {
                logger.warn("Input bam was marked as " + mergedHeader.getSortOrder().toString() + " but " + TREAT_UNSORTED_AS_ORDERED + " is specified so it's being treated as read name grouped");
                mergedHeader.setGroupOrder(SAMFileHeader.GroupOrder.query);
            }
        }

        JavaRDD<GATKRead> reads = getReads();
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null);
        // If we need to remove optical duplicates, set the engine to mark optical duplicates using the DT tag.
        if (markDuplicatesSparkArgumentCollection.removeSequencingDuplicates && markDuplicatesSparkArgumentCollection.taggingPolicy == MarkDuplicates.DuplicateTaggingPolicy.DontTag) {
            markDuplicatesSparkArgumentCollection.taggingPolicy = MarkDuplicates.DuplicateTaggingPolicy.OpticalOnly;
        }

        final JavaRDD<GATKRead> finalReadsForMetrics = mark(reads, mergedHeader, finder, markDuplicatesSparkArgumentCollection, getRecommendedNumReducers());

        if (metricsFile != null) {
            final JavaPairRDD<String, GATKDuplicationMetrics> metricsByLibrary = MarkDuplicatesSparkUtils.generateMetrics(
                    mergedHeader, finalReadsForMetrics);
            final MetricsFile<GATKDuplicationMetrics, Double> resultMetrics = getMetricsFile();
            MarkDuplicatesSparkUtils.saveMetricsRDD(resultMetrics, mergedHeader, metricsByLibrary, metricsFile);
        }
        JavaRDD<GATKRead> readsForWriting = finalReadsForMetrics;
        // Filter out the duplicates if instructed to do so
        if (markDuplicatesSparkArgumentCollection.removeAllDuplicates) {
            readsForWriting = readsForWriting.filter(r -> !r.isDuplicate());
        } else if (markDuplicatesSparkArgumentCollection.removeSequencingDuplicates) {
            readsForWriting = readsForWriting.filter(r -> !MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(r.getAttributeAsString(MarkDuplicates.DUPLICATE_TYPE_TAG)));
        }

        mergedHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        writeReads(ctx, output, readsForWriting, mergedHeader, true);
    }

    // helper method to determin if an input header is to be treated as a query group sorted file.
    private boolean treatAsReadGroupOrdered(SAMFileHeader header, boolean treatUnsortedAsReadGrouped) {
        final SAMFileHeader.SortOrder sortOrder = header.getSortOrder();
        if( ReadUtils.isReadNameGroupedBam(header) ){
            return true;
        } else if ( treatUnsortedAsReadGrouped && (sortOrder.equals(SAMFileHeader.SortOrder.unknown) || sortOrder.equals(SAMFileHeader.SortOrder.unsorted))) {
            logger.warn("Input bam was marked as " + sortOrder.toString() + " but " + TREAT_UNSORTED_AS_ORDERED + " is specified so it's being treated as read name grouped");
            return true;
        }
        return false;
    }

}
