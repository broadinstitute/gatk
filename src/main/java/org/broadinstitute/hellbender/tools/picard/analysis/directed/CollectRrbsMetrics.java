package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.util.*;

/**
 * Calculates and reports QC metrics for RRBS data based on the methylation status at individual C/G bases as well
 * as CpG sites across all reads in the input BAM/SAM file.
 *
 * @author jgentry@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = CollectRrbsMetrics.USAGE,
        usageShort = CollectRrbsMetrics.USAGE,
        programGroup = QCProgramGroup.class
)
public final class CollectRrbsMetrics extends PicardCommandLineProgram {
    final static String USAGE = "Collects metrics about bisulfite conversion for RRBS data";
    private static final String R_SCRIPT = "rrbsQc.R";

    @Argument(doc = "The BAM or SAM file containing aligned reads. Must be coordinate sorted", shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Base name for output files", shortName = StandardArgumentDefinitions.METRICS_FILE_SHORT_NAME)
    public String METRICS_FILE_PREFIX;

    @Argument(doc = "Minimum read length")
    public int MINIMUM_READ_LENGTH = 5;

    @Argument(doc = "Threshold for base quality of a C base before it is considered")
    public int C_QUALITY_THRESHOLD = 20;

    @Argument(doc = "Threshold for quality of a base next to a C before the C base is considered")
    public int NEXT_BASE_QUALITY_THRESHOLD = 10;

    @Argument(doc = "Maximum percentage of mismatches in a read for it to be considered, with a range of 0-1")
    public double MAX_MISMATCH_RATE = 0.1;

    @Argument(doc = "Set of sequence names to consider, if not specified all sequences will be used", optional = true)
    public Set<String> SEQUENCE_NAMES = new HashSet<>();

    @Argument(shortName = StandardArgumentDefinitions.ASSUME_SORTED_SHORT_NAME,
            doc = "If true, assume that the input file is coordinate sorted even if the header says otherwise.")
    public boolean ASSUME_SORTED = false;

    @Argument(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Argument(doc = "Should an output plot be created")
    public boolean PRODUCE_PLOT = false;

    public static final String DETAIL_FILE_EXTENSION = "rrbs_detail_metrics";
    public static final String SUMMARY_FILE_EXTENSION = "rrbs_summary_metrics";
    public static final String PDF_FILE_EXTENSION = "rrbs_qc.pdf";

    private static final Log log = Log.getInstance(CollectRrbsMetrics.class);

    @Override
    protected Object doWork() {
        if (!METRICS_FILE_PREFIX.endsWith(".")) {
            METRICS_FILE_PREFIX = METRICS_FILE_PREFIX + ".";
        }
        final File SUMMARY_OUT = new File(METRICS_FILE_PREFIX + SUMMARY_FILE_EXTENSION);
        final File DETAILS_OUT = new File(METRICS_FILE_PREFIX + DETAIL_FILE_EXTENSION);
        final File PLOTS_OUT = new File(METRICS_FILE_PREFIX + PDF_FILE_EXTENSION);
        assertIoFiles(SUMMARY_OUT, DETAILS_OUT, PLOTS_OUT);

        final SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        if (!ASSUME_SORTED && samReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new UserException("The input file " + INPUT.getAbsolutePath() + " does not appear to be coordinate sorted");
        }

        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final ProgressLogger progressLogger = new ProgressLogger(log);

        final RrbsMetricsCollector metricsCollector = new RrbsMetricsCollector(METRIC_ACCUMULATION_LEVEL, samReader.getFileHeader().getReadGroups(),
                C_QUALITY_THRESHOLD, NEXT_BASE_QUALITY_THRESHOLD, MINIMUM_READ_LENGTH, MAX_MISMATCH_RATE);

        for (final SAMRecord samRecord : samReader) {
            progressLogger.record(samRecord);
            if (!samRecord.getReadUnmappedFlag() && !isSequenceFiltered(samRecord.getReferenceName())) {
                final ReferenceSequence referenceSequence = refWalker.get(samRecord.getReferenceIndex());
                metricsCollector.acceptRecord(samRecord, referenceSequence);
            }
        }
        metricsCollector.finish();
        final MetricsFile<RrbsMetrics, Long> rrbsMetrics = getMetricsFile();
        metricsCollector.addAllLevelsToFile(rrbsMetrics);

        // Using RrbsMetrics as a way to get both of the metrics objects through the MultiLevelCollector. Once
        // we get it out split it apart to the two separate MetricsFiles and write them to file
        final MetricsFile<RrbsSummaryMetrics, ?> summaryFile = getMetricsFile();
        final MetricsFile<RrbsCpgDetailMetrics, ?> detailsFile = getMetricsFile();
        for (final RrbsMetrics rrbsMetric : rrbsMetrics.getMetrics()) {
            summaryFile.addMetric(rrbsMetric.getSummaryMetrics());
            for (final RrbsCpgDetailMetrics detailMetric : rrbsMetric.getDetailMetrics()) {
                detailsFile.addMetric(detailMetric);
            }
        }
        summaryFile.write(SUMMARY_OUT);
        detailsFile.write(DETAILS_OUT);
        if(PRODUCE_PLOT) {
            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(new Resource(R_SCRIPT, CollectRrbsMetrics.class));
            executor.addArgs(DETAILS_OUT.getAbsolutePath(), SUMMARY_OUT.getAbsolutePath(), PLOTS_OUT.getAbsolutePath());
            executor.exec();
        }

        CloserUtil.close(samReader);
        return null;
    }

    private boolean isSequenceFiltered(final String sequenceName) {
        return (SEQUENCE_NAMES != null) && (SEQUENCE_NAMES.size() > 0) && (!SEQUENCE_NAMES.contains(sequenceName));
    }

    private void assertIoFiles(final File summaryFile, final File detailsFile, final File plotsFile) {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsWritable(summaryFile);
        IOUtil.assertFileIsWritable(detailsFile);
        IOUtil.assertFileIsWritable(plotsFile);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errorMsgs = new ArrayList<>();
        if (MAX_MISMATCH_RATE < 0 || MAX_MISMATCH_RATE > 1) {
            errorMsgs.add("MAX_MISMATCH_RATE must be in the range of 0-1");
        }

        if (C_QUALITY_THRESHOLD < 0) {
            errorMsgs.add("C_QUALITY_THRESHOLD must be >= 0");
        }

        if (NEXT_BASE_QUALITY_THRESHOLD < 0) {
            errorMsgs.add("NEXT_BASE_QUALITY_THRESHOLD must be >= 0");
        }

        if (MINIMUM_READ_LENGTH <= 0) {
            errorMsgs.add("MINIMUM_READ_LENGTH must be > 0");
        }

        return errorMsgs.size() == 0 ? null : errorMsgs.toArray(new String[errorMsgs.size()]);
    }
}
