package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.metrics.*;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.File;

/**
 * Command line program to read non-duplicate insert sizes, create a Histogram
 * and report distribution statistics.
 */
@CommandLineProgramProperties(
        summary = "Reads a SAM/BAM/CRAM file and writes a file containing metrics about " +
                "the statistical distribution of insert size (excluding duplicates) " +
                "and generates a Histogram plot.",
        oneLineSummary = "Produces metrics for insert size distribution for a SAM/BAM/CRAM file",
        programGroup = QCProgramGroup.class
)
public final class CollectInsertSizeMetrics extends SinglePassSamProgram {
    private static final Logger log = LogManager.getLogger(CollectInsertSizeMetrics.class);

    @ArgumentCollection
    public InsertSizeMetricsArgumentCollection insertSizeArgs = new InsertSizeMetricsArgumentCollection();

    // Calculates InsertSizeMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private InsertSizeMetricsCollector insertSizeCollector;

    private ReadFilter insertSizeMetricsReadFilter = null;

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
         if (insertSizeArgs.minimumPct < 0 || insertSizeArgs.minimumPct > 0.5) {
             return new String[]{"MINIMUM_PCT was set to " + insertSizeArgs.minimumPct +
                     ". It must be between 0 and 0.5 so all data categories don't get discarded."};
         }

         return super.customCommandLineValidation();
    }

    @Override
    protected boolean usesNoRefReads() { return false; }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(new File(insertSizeArgs.output));
        if (insertSizeArgs.producePlot) {
            IOUtil.assertFileIsWritable(new File(insertSizeArgs.histogramPlotFile));
        }

        //Delegate actual collection to InsertSizeMetricCollector
        insertSizeCollector = new InsertSizeMetricsCollector(insertSizeArgs, header);
        insertSizeMetricsReadFilter = insertSizeCollector.getReadFilter(header);
    }

    @Override
    protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
        // InsertSizeMetricsCollector assumes that any records passed have already been filtered
        if (insertSizeMetricsReadFilter.test(new SAMRecordToGATKReadAdapter(record))) {
            insertSizeCollector.acceptRecord(record, ref);
        }
    }

    @Override
    protected void finish() {
        final MetricsFile<InsertSizeMetrics, Integer> metricsFile = getMetricsFile();
        insertSizeCollector.finish(metricsFile, INPUT.getName(), null);
    }
}
