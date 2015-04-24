package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;

import java.io.File;

/**
 * Command line program to calibrate quality yield metrics
 *
 * @author Martha Borkan
 */
@CommandLineProgramProperties(
        usage = "Collects quality yield metrics, a set of metrics that quantify the quality and yield of sequence data from a " +
                "SAM/BAM input file.",
        usageShort = "Collects a set of metrics that quantify the quality and yield of sequence data from the provided SAM/BAM",
        programGroup = QCProgramGroup.class
)
public final class CollectQualityYieldMetrics extends PicardCommandLineProgram {

    @Argument(shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "A SAM or BAM file to process.")
    public File INPUT;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The metrics file to write with quality yield metrics.")
    public File OUTPUT;

    @Argument(shortName = StandardArgumentDefinitions.USE_ORIGINAL_QUALITIES_SHORT_NAME,
            doc = "If available in the OQ tag, use the original quality scores " +
                    "as inputs instead of the quality scores in the QUAL field.")
    public boolean USE_ORIGINAL_QUALITIES = true;

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then iterates through
     * all the records accumulating metrics.  Finally writes metrics file
     */
    protected Object doWork() {
        final Log log = Log.getInstance(getClass());
        final ProgressLogger progress = new ProgressLogger(log);

        // Some quick parameter checking
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        log.info("Reading input file and calculating metrics.");

        final SamReader sam = SamReaderFactory.makeDefault().open(INPUT);

        final MetricsFile<QualityYieldMetrics, Integer> metricsFile = getMetricsFile();
        final QualityYieldMetrics metrics = new QualityYieldMetrics();

        for (final SAMRecord rec : sam) {
            metrics.TOTAL_READS++;
            final int length = rec.getReadLength();

            final boolean isPfRead = !rec.getReadFailsVendorQualityCheckFlag();
            if (isPfRead) {
                metrics.PF_READS++;
                metrics.PF_BASES += length;
            }

            metrics.TOTAL_BASES += length;

            final byte[] quals;
            if (USE_ORIGINAL_QUALITIES) {
                byte[] tmp = rec.getOriginalBaseQualities();
                if (tmp == null) tmp = rec.getBaseQualities();
                quals = tmp;
            } else {
                quals = rec.getBaseQualities();
            }

            // add up quals, and quals >= 20
            for (int i = 0; i < quals.length; ++i) {
                metrics.Q20_EQUIVALENT_YIELD += quals[i];
                if (quals[i] >= 20) metrics.Q20_BASES++;
                if (quals[i] >= 30) metrics.Q30_BASES++;

                if (isPfRead) {
                    metrics.PF_Q20_EQUIVALENT_YIELD += quals[i];
                    if (quals[i] >= 20) metrics.PF_Q20_BASES++;
                    if (quals[i] >= 30) metrics.PF_Q30_BASES++;
                }
            }

            progress.record(rec);
        }

        metrics.READ_LENGTH = metrics.TOTAL_READS == 0 ? 0 : (int) (metrics.TOTAL_BASES / metrics.TOTAL_READS);
        metrics.Q20_EQUIVALENT_YIELD = metrics.Q20_EQUIVALENT_YIELD / 20;
        metrics.PF_Q20_EQUIVALENT_YIELD = metrics.PF_Q20_EQUIVALENT_YIELD / 20;

        metricsFile.addMetric(metrics);
        metricsFile.write(OUTPUT);

        return null;
    }

    /** A set of metrics used to describe the general quality of a BAM file */
    public static class QualityYieldMetrics extends MetricBase {

        /** The total number of reads in the input file */
        public int TOTAL_READS = 0;

        /** The number of reads that are PF - pass filter */
        public int PF_READS = 0;

        /** The average read length of all the reads (will be fixed for a lane) */
        public int READ_LENGTH = 0;

        /** The total number of bases in all reads */
        public long TOTAL_BASES;

        /** The total number of bases in all PF reads */
        public long PF_BASES = 0;

        /** The number of bases in all reads that achieve quality score 20 or higher */
        public long Q20_BASES = 0;

        /** The number of bases in PF reads that achieve quality score 20 or higher */
        public long PF_Q20_BASES = 0;

        /** The number of bases in all reads that achieve quality score 20 or higher */
        public long Q30_BASES = 0;

        /** The number of bases in PF reads that achieve quality score 20 or higher */
        public long PF_Q30_BASES = 0;

        /** The sum of quality scores of all bases divided by 20 */
        public long Q20_EQUIVALENT_YIELD = 0;

        /** The sum of quality scores of all bases divided by 20 */
        public long PF_Q20_EQUIVALENT_YIELD = 0;

    }

}
