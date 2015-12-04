package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;

/**
 * Command line program to calibrate quality yield metrics
 *
 * @author Martha Borkan
 */
@CommandLineProgramProperties(
        summary = "Collects quality yield metrics, a set of metrics that quantify the quality and yield of sequence data from a " +
                "SAM/BAM/CRAM input file.",
        oneLineSummary = "Produces metrics that quantify the quality and yield of sequence data from the provided SAM/BAM/CRAM file",
        programGroup = QCProgramGroup.class
)
public final class CollectQualityYieldMetrics extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "A SAM/BAM/CRAM file to process.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The metrics file to write with quality yield metrics.")
    public File OUTPUT;

    @Argument(fullName = StandardArgumentDefinitions.USE_ORIGINAL_QUALITIES_LONG_NAME, shortName = StandardArgumentDefinitions.USE_ORIGINAL_QUALITIES_SHORT_NAME,
            doc = "If available in the OQ tag, use the original quality scores " +
                    "as inputs instead of the quality scores in the QUAL field.")
    public boolean USE_ORIGINAL_QUALITIES = true;

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then iterates through
     * all the records accumulating metrics.  Finally writes metrics file
     */
    protected Object doWork() {
        final ProgressLogger progress = new ProgressLogger(logger);

        // Some quick parameter checking
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        logger.info("Reading input file and calculating metrics.");

        final SamReader sam = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

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

        /** The number of bases in all reads that achieve quality score 30 or higher */
        public long Q30_BASES = 0;

        /** The number of bases in PF reads that achieve quality score 30 or higher */
        public long PF_Q30_BASES = 0;

        /** The sum of quality scores of all bases divided by 20 */
        public long Q20_EQUIVALENT_YIELD = 0;

        /** The sum of quality scores of all bases divided by 20 */
        public long PF_Q20_EQUIVALENT_YIELD = 0;

    }

}
