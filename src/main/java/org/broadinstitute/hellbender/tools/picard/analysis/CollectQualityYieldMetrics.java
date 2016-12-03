package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.metrics.QualityYieldMetrics;
import org.broadinstitute.hellbender.metrics.QualityYieldMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
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

    @ArgumentCollection
    public QualityYieldMetricsArgumentCollection qualityYieldMetricsArgs = new QualityYieldMetricsArgumentCollection();

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then iterates through
     * all the records accumulating metrics.  Finally writes metrics file
     */
    @Override
    protected Object doWork() {
        final ProgressLogger progress = new ProgressLogger(logger);

        File output = new File(qualityYieldMetricsArgs.output);
        // Some quick parameter checking
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(output);

        logger.info("Reading input file and calculating metrics.");

        final SamReader sam = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final MetricsFile<QualityYieldMetrics, Integer> metricsFile = getMetricsFile();
        final QualityYieldMetrics metrics = new QualityYieldMetrics();

        metrics.setUseOriginalQualities(qualityYieldMetricsArgs.useOriginalQualities);

        for (final SAMRecord rec : sam) {
            metrics.addRead(new SAMRecordToGATKReadAdapter(rec));
            progress.record(rec);
        }

        metrics.finish();
        metricsFile.addMetric(metrics);
        metricsFile.write(output);

        return null;
    }

}
