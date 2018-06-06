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
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.GATKDuplicationMetrics;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.markduplicates.MarkDuplicates;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import scala.Tuple2;

import java.util.*;

/**
 * <p>This is a Spark implementation of the MarkDuplicates tool from Picard that allows the tool to be run in
 *    parallel on multiple cores on a local machine or multiple machines on a Spark cluster while still matching
 *    the output of the single-core Picard version. Since the tool requires holding all of the readnames in memory
 *    while it groups the read information, it is recommended running this tool on a machine/configuration
 *    with at least 8 GB of memory overall for a typical 30x bam.</p>
 *
 * <p>This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are
 *    defined as originating from a single fragment of DNA.  Duplicates can arise during sample preparation e.g. library
 *    construction using PCR.  See also "<a href='https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity'>EstimateLibraryComplexity</a>" +
 *    for additional notes on PCR duplication artifacts.  Duplicate reads can also result from a single amplification cluster,
 *    incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.  These duplication artifacts are
 *    referred to as optical duplicates.</p>
 *
 * <p>The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file.
 *    After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks
 *    reads by the sums of their base-quality scores (default method).</p>
 *
 * <p>The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each
 *    read.  Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024.
 *    If you are not familiar with this type of annotation, please see the following <a href='https://www.broadinstitute.org/gatk/blog?id=7019'>blog post</a> for additional information.</p>" +
 *
 * <p>Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of
 *    duplicate.  To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in
 *    the 'optional field' section of a SAM/BAM file.  Invoking the 'duplicate-tagging-policy' option,
 *    you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no
 *    duplicates (DontTag).  The records within the output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked
 *    'duplicate-tagging-policy'), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ).
 *    This tool uses the 'read-name-regex' and the 'optical-duplicate-pixel-distance' options as the primary methods to identify
 *    and differentiate duplicate types.  Set read-name-regex' to null to skip optical duplicate detection, e.g. for RNA-seq
 *    or other data where duplicate sets are extremely large and estimating library complexity is not an aim.
 *    Note that without optical duplicate counts, library size estimation will be inaccurate.</p>
 *
 * <p>MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.</p>
 *
 * <p>The program can take either coordinate-sorted or query-sorted inputs, however it is recommended that the input be
 *    query-sorted or query-grouped as the tool will have to perform an extra sort operation on the data in order to associate
 *    reads from the input bam with their mates.</p>
 *
 * <p>If desired, duplicates can be removed using the 'remove-all-duplicates' and 'remove-sequencing-duplicates' options.</p>
 *
 * <h4>Usage example:</h4>
 *     <pre>
 *      gatk MarkDuplicatesSpark \\<br />
 *            -I input.bam \\<br />
 *            -O marked_duplicates.bam \\<br />
 *            -M marked_dup_metrics.txt
 *     </pre>
 *
 *  <h4>MarkDuplicates run locally specifying the core input (if 'spark.executor.cores' is unset spark will use all available cores on the machine)</h4>
 *     <pre>
 *       gatk MarkDuplicatesSpark \\<br />
 *            -I input.bam \\<br />
 *            -O marked_duplicates.bam \\<br />
 *            -M marked_dup_metrics.txt \\<br />
 *            --conf 'spark.executor.cores=5'
 *     </pre>
 *
 *  <h4>MarkDuplicates run on a spark cluster 5 machines</h4>
 *     <pre>
 *       gatk MarkDuplicatesSpark \\<br />
 *            -I input.bam \\<br />
 *            -O marked_duplicates.bam \\<br />
 *            -M marked_dup_metrics.txt \\<br />
 *            -- \\<br />
 *            --spark-runner SPARK \\<br />
 *            --spark-master <master_url> \\<br />
 *            --num-executors 5 \\<br />
 *            --executor-cores 8 <br />
 *     </pre>
 *
 *    Please see
 *    <a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics'>MarkDuplicates</a>
 *    for detailed explanations of the output metrics.
 *    <hr />
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary ="Marks duplicates on spark",
        oneLineSummary ="MarkDuplicates on Spark",
        programGroup = ReadDataManipulationProgramGroup.class)
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
                                        ((SAMRecordToGATKReadAdapter) read).setTransientAttribute(MarkDuplicatesSparkUtils.OPTICAL_DUPLICATE_TOTAL_ATTRIBUTE_NAME, dupCount);
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
        // Check if we are using multiple inputs that the headers are all in the correct querygrouped ordering
        Map<String, SAMFileHeader> headerMap = getReadSouceHeaderMap();
        if (headerMap.size() > 1) {
            headerMap.entrySet().stream().forEach(h -> {if(!ReadUtils.isReadNameGroupedBam(h.getValue())) {
                throw new UserException("Multiple inputs to MarkDuplicatesSpark detected but input "+h.getKey()+" was sorted in "+h.getValue().getSortOrder()+" order");
                    }});
        }

        JavaRDD<GATKRead> reads = getReads();
        final OpticalDuplicateFinder finder = opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null ?
                new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;
        // If we need to remove optical duplicates, set the engine to mark optical duplicates using the DT tag.
        if (markDuplicatesSparkArgumentCollection.removeSequencingDuplicates && markDuplicatesSparkArgumentCollection.taggingPolicy == MarkDuplicates.DuplicateTaggingPolicy.DontTag) {
            markDuplicatesSparkArgumentCollection.taggingPolicy = MarkDuplicates.DuplicateTaggingPolicy.OpticalOnly;
        }

        final SAMFileHeader header = getHeaderForReads();
        final JavaRDD<GATKRead> finalReadsForMetrics = mark(reads, header, finder, markDuplicatesSparkArgumentCollection, getRecommendedNumReducers());

        if (metricsFile != null) {
            final JavaPairRDD<String, GATKDuplicationMetrics> metricsByLibrary = MarkDuplicatesSparkUtils.generateMetrics(
                    header, finalReadsForMetrics);
            final MetricsFile<GATKDuplicationMetrics, Double> resultMetrics = getMetricsFile();
            MarkDuplicatesSparkUtils.saveMetricsRDD(resultMetrics, header, metricsByLibrary, metricsFile);
        }
        JavaRDD<GATKRead> readsForWriting = finalReadsForMetrics;
        // Filter out the duplicates if instructed to do so
        if (markDuplicatesSparkArgumentCollection.removeAllDuplicates) {
            readsForWriting = readsForWriting.filter(r -> !r.isDuplicate());
        } else if (markDuplicatesSparkArgumentCollection.removeSequencingDuplicates) {
            readsForWriting = readsForWriting.filter(r -> !MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(r.getAttributeAsString(MarkDuplicates.DUPLICATE_TYPE_TAG)));
        }

        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        writeReads(ctx, output, readsForWriting, header, true);
    }

}
