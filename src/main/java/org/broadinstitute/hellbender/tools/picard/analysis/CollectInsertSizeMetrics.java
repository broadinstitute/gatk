package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.util.EnumSet;
import java.util.Set;

/**
 * Command line program to read non-duplicate insert sizes, create a Histogram
 * and report distribution statistics.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
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

    private static final String R_SCRIPT = "insertSizeHistogram.R";

    @Argument(shortName="H", doc="File to write insert size Histogram chart to.")
    public File HISTOGRAM_FILE;

    @Argument(doc="Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. " +
            "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
            "artifacts to make the mean and sd grossly misleading regarding the real distribution.")
    public double DEVIATIONS = 10;

    @Argument(shortName="W", doc="Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. " +
            "Also, when calculating mean and standard deviation, only bins <= HISTOGRAM_WIDTH will be included.", optional=true)
    public Integer HISTOGRAM_WIDTH = null;

    @Argument(shortName="M", doc="When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this " +
            "percentage of overall reads. (Range: 0 to 1).")
    public float MINIMUM_PCT = 0.05f;

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = EnumSet.of(MetricAccumulationLevel.ALL_READS);

    @Argument(doc = "Should an output plot be created")
    public boolean PRODUCE_PLOT = false;

    // Calculates InsertSizeMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private InsertSizeMetricsCollector multiCollector;

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
         if (MINIMUM_PCT < 0 || MINIMUM_PCT > 0.5) {
             return new String[]{"MINIMUM_PCT was set to " + MINIMUM_PCT + ". It must be between 0 and 0.5 so all data categories don't get discarded."};
         }

         return super.customCommandLineValidation();
    }

    @Override
    protected boolean usesNoRefReads() { return false; }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(HISTOGRAM_FILE);

        //Delegate actual collection to InsertSizeMetricCollector
        multiCollector = new InsertSizeMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), MINIMUM_PCT, HISTOGRAM_WIDTH, DEVIATIONS);
    }

    @Override
    protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
        multiCollector.acceptRecord(record, ref);
    }

    @Override
    protected void finish() {
        multiCollector.finish();

        final MetricsFile<InsertSizeMetrics, Integer> file = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);

        if(file.getNumHistograms() == 0) {
            //can happen if user sets MINIMUM_PCT = 0.5, etc.
            log.warn("All data categories were discarded because they contained < " + MINIMUM_PCT +
                     " of the total aligned paired data.");
            final InsertSizeMetricsCollector.PerUnitInsertSizeMetricsCollector allReadsCollector
                    = (InsertSizeMetricsCollector.PerUnitInsertSizeMetricsCollector) multiCollector.getAllReadsCollector();
            log.warn("Total mapped pairs in all categories: " + (allReadsCollector == null ? allReadsCollector : allReadsCollector.getTotalInserts()));
        }
        else  {
            file.write(OUTPUT);
            if(PRODUCE_PLOT){
                final RScriptExecutor executor = new RScriptExecutor();
                executor.addScript(new Resource(R_SCRIPT, CollectInsertSizeMetrics.class));
                executor.addArgs(OUTPUT.getAbsolutePath(), HISTOGRAM_FILE.getAbsolutePath(), INPUT.getName());
                if (HISTOGRAM_WIDTH != null) executor.addArgs(String.valueOf(HISTOGRAM_WIDTH));
                executor.exec();
            }
        }
    }
}
