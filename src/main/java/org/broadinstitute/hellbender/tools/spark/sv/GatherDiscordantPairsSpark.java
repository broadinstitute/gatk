package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import scala.Tuple2;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

@CommandLineProgramProperties(summary="Gather discordant read pairs using spark",
        oneLineSummary="Gather discordant read pairs using spark",
        programGroup = SparkProgramGroup.class)
public class GatherDiscordantPairsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Argument(doc = "picard insert size metrics file", shortName = "picardIsMetricsFile",
            fullName = "picardIsMetricsFile", optional = false)
    private String picardInsertSizeMetricsFile;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        try {
            final InsertSizeDistribution insertSizeDistribution = parsePicardInsertSizeMetricsFile(picardInsertSizeMetricsFile);

            final JavaRDD<GATKRead> reads = getReads();   // get all reads passing WellFormedReadFilter
            final SAMFileHeader header = getHeaderForReads();
            final JavaRDD<GATKRead> discordantReadPairs = DiscordantPairFilterFn.getDiscordantReadPairs(insertSizeDistribution, reads, header);
            try {
                ReadsSparkSink.writeReads(ctx,output,discordantReadPairs,header, ReadsWriteFormat.SINGLE); }
            catch ( IOException e ) {
                throw new GATKException("Unable to write BAM" + output, e);
            }

        } catch (FileNotFoundException e) {
            throw new UserException("Unable to read insert size metrics file " + picardInsertSizeMetricsFile, e);
        }

    }

    private InsertSizeDistribution parsePicardInsertSizeMetricsFile(final String picardInsertSizeMetricsFile) throws FileNotFoundException {
        MetricsFile<InsertSizeMetrics, Long> metricsFile = new MetricsFile<InsertSizeMetrics, Long>();
        Reader picardInsertSizeMetricsReader = new BufferedReader(new FileReader(picardInsertSizeMetricsFile));
        metricsFile.read(picardInsertSizeMetricsReader);

        return new PicardInsertSizeMetrics(metricsFile);
    }

    private static interface InsertSizeDistribution {
        public int getMedian(String sample, String library, String readGroup);

        public int getMAD(String sample, String library, String readGroup);

        boolean hasInfoFor(String sample, String library, String id);
    }

    /**
     * Read metrics from a metrics file that was produced by Picard::CollectInsertSizeMetrics,
     *   results are grouped by {SAMPLE, LIBRARY, READ_GROUP}
     */
    private static class PicardInsertSizeMetrics implements InsertSizeDistribution, Serializable {
        private static final long serialVersionUID = 1L;

        Map<String, Double> medianBySample = new HashMap<>();
        Map<String, Double> madBySample = new HashMap<>();
        Map<String, Double> medianByLibrary = new HashMap<>();
        Map<String, Double> madByLibrary = new HashMap<>();
        Map<String, Double> medianByReadGroup = new HashMap<>();
        Map<String, Double> madByReadGroup = new HashMap<>();

        public PicardInsertSizeMetrics(final MetricsFile<InsertSizeMetrics, Long> picardMetricsFile) {
            for (int i = 0; i < picardMetricsFile.getMetrics().size(); i++) {
                InsertSizeMetrics insertSizeMetrics = picardMetricsFile.getMetrics().get(i);
                if (insertSizeMetrics.READ_GROUP != null) {
                    medianByReadGroup.put(insertSizeMetrics.READ_GROUP, insertSizeMetrics.MEDIAN_INSERT_SIZE);
                    madByReadGroup.put(insertSizeMetrics.READ_GROUP, insertSizeMetrics.MEDIAN_ABSOLUTE_DEVIATION);
                    continue;
                }
                if (insertSizeMetrics.LIBRARY != null) {
                    medianByLibrary.put(insertSizeMetrics.LIBRARY, insertSizeMetrics.MEDIAN_INSERT_SIZE);
                    madByLibrary.put(insertSizeMetrics.LIBRARY, insertSizeMetrics.MEDIAN_ABSOLUTE_DEVIATION);
                    continue;
                }
                if (insertSizeMetrics.SAMPLE != null) {
                    medianBySample.put(insertSizeMetrics.SAMPLE, insertSizeMetrics.MEDIAN_INSERT_SIZE);
                    madBySample.put(insertSizeMetrics.SAMPLE, insertSizeMetrics.MEDIAN_ABSOLUTE_DEVIATION);
                    continue;
                }
            }
        }

        @Override
        public int getMedian(final String sample, final String library, final String readGroup) {
            if (medianByReadGroup.containsKey(readGroup)) {
                return medianByReadGroup.get(readGroup).intValue();
            }
            if (medianByLibrary.containsKey(library)) {
                return medianByLibrary.get(library).intValue();
            }
            if (medianBySample.containsKey(sample)) {
                return medianBySample.get(sample).intValue();
            }
            throw new UserException("No insert size metrics found for sample/library/readgroup " + sample + "/" + library + "/" + readGroup);
        }

        @Override
        public int getMAD(final String sample, final String library, final String readGroup) {
            if (madByReadGroup.containsKey(readGroup)) {
                return madByReadGroup.get(readGroup).intValue();
            }
            if (madByLibrary.containsKey(library)) {
                return madByLibrary.get(library).intValue();
            }
            if (madBySample.containsKey(sample)) {
                return madBySample.get(sample).intValue();
            }
            throw new UserException("No insert size metrics found for sample/library/readgroup " + sample + "/" + library + "/" + readGroup);
        }

        @Override
        public boolean hasInfoFor(final String sample, final String library, final String readGroup) {
            return medianBySample.containsKey(sample) || medianByLibrary.containsKey(library) || medianByReadGroup.containsKey(readGroup);
        }
    }

    private static class DiscordantPairFilterFn implements Serializable {
        private static final long serialVersionUID = 1L;
        private final SAMFileHeader header;
        private final InsertSizeDistribution insertSizeDistribution;

        /**
         * Gather discordant read pairs based on "null distribution" information gathered from Picard-gernated metrics file.
         * Three cases are tested: case 1: pairs not oriented FR
         *                         case 2: ends map to different chromosome
         *                         case 3: out of normal range sized pairs
         * @param insertSizeDistribution  "null distribution"
         * @param reads                   reads to be tested
         * @param header                  header in input file
         * @return                        discordant paired reads
         */
        public static JavaRDD<GATKRead> getDiscordantReadPairs(final InsertSizeDistribution insertSizeDistribution, final JavaRDD<GATKRead> reads, final SAMFileHeader header) {
            return reads
                    .filter(read ->
                            new DiscordantPairFilterFn(header, insertSizeDistribution).isPartOfDiscordantPair(read))
                    .mapToPair(read -> new Tuple2<>(read.getName(), read))
                    .groupByKey()
                    .flatMap(Tuple2::_2);
        }

        public DiscordantPairFilterFn(final SAMFileHeader header, final InsertSizeDistribution insertSizeDistribution) {
            this.header = header;
            this.insertSizeDistribution = insertSizeDistribution;
        }

        private boolean isPartOfDiscordantPair(GATKRead read) {

            // filter reads not to be considered
            if (read.failsVendorQualityCheck() ||
                    read.isDuplicate() ||
                    read.isSecondaryAlignment() ||
                    read.isSupplementaryAlignment() ||
                    read.mateIsUnmapped() ||
                    read.getMappingQuality() == 0 ||
                    read.getAttributeAsInteger("MQ") == 0) {
                return false;
            }

            final SAMReadGroupRecord readGroup = header.getReadGroup(read.getReadGroup());

            if (!insertSizeDistribution.hasInfoFor(readGroup.getSample(), readGroup.getLibrary(), readGroup.getId())) {
                return false;
            }

            // case 1: pairs not oriented FR
            if (! (SamPairUtil.getPairOrientation(read.convertToSAMRecord(header)) == SamPairUtil.PairOrientation.FR)) {
                return true;
            }

            // case 2: ends map to different chromosome
            if (! read.getContig().equals(read.getMateContig())) {
                return true;
            }

            // case 3: out of normal range sized pairs
            final int isizeThreshold = (insertSizeDistribution.getMedian(readGroup.getSample(), readGroup.getLibrary(), readGroup.getId())) +
                    3 * insertSizeDistribution.getMAD(readGroup.getSample(), readGroup.getLibrary(), readGroup.getId());

            if (Math.abs(read.getFragmentLength()) > isizeThreshold) {
                return true;
            }

            return false;
        }
    }
}
