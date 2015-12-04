package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.util.*;

/**
 * Program to generate a data table and chart of mean quality by cycle from a
 * BAM file.  Works best on a single lane/run of data, but can be applied to
 * merged BAMs - the output may just be a little confusing.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Program to generate a data table and pdf chart of " +
                "mean base quality by cycle from a SAM/BAM/CRAM file.  Works best on a single lane/run of data, but can be applied to" +
                "merged BAMs. Uses R to generate chart output.",
        oneLineSummary = "Produces metrics for mean quality by cycle for a SAM/BAM/CRAM file",
        programGroup = QCProgramGroup.class
)
public final class MeanQualityByCycle extends SinglePassSamProgram {
    public static final String R_SCRIPT = "meanQualityByCycle.R";

    @Argument(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.", optional = true)
    public File CHART_OUTPUT;

    @Argument(doc="If set to true, calculate mean quality over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(doc="If set to true calculate mean quality over PF reads only.")
    public boolean PF_READS_ONLY = false;

    @Argument(doc = "Should an output plot be created")
    public boolean PRODUCE_PLOT = false;

    private final HistogramGenerator q  = new HistogramGenerator(false);
    private final HistogramGenerator oq = new HistogramGenerator(true);

    /**
     * A subtitle for the plot, usually corresponding to a library.
     */
    private String plotSubtitle = "";

    private final Logger logger = LogManager.getLogger(this.getClass());

    private static class HistogramGenerator {
        final boolean useOriginalQualities;
        int maxLengthSoFar = 0;
        double[] firstReadTotalsByCycle  = new double[maxLengthSoFar];
        long[]   firstReadCountsByCycle  = new long[maxLengthSoFar];
        double[] secondReadTotalsByCycle = new double[maxLengthSoFar];
        long[]   secondReadCountsByCycle = new long[maxLengthSoFar];

        private HistogramGenerator(final boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        void addRecord(final SAMRecord rec) {
            final byte[] quals = (useOriginalQualities ? rec.getOriginalBaseQualities() : rec.getBaseQualities());
            if (quals == null) return;

            final int length = quals.length;
            final boolean rc = rec.getReadNegativeStrandFlag();
            ensureArraysBigEnough(length+1);

            for (int i=0; i<length; ++i) {
                final int cycle = rc ? length-i : i+1;

                if (rec.getReadPairedFlag() && rec.getSecondOfPairFlag()) {
                    secondReadTotalsByCycle[cycle] += quals[i];
                    secondReadCountsByCycle[cycle] += 1;
                }
                else {
                    firstReadTotalsByCycle[cycle] += quals[i];
                    firstReadCountsByCycle[cycle] += 1;
                }
            }
        }

        private void ensureArraysBigEnough(final int length) {
            if (length > maxLengthSoFar) {
                firstReadTotalsByCycle  = Arrays.copyOf(firstReadTotalsByCycle, length);
                firstReadCountsByCycle  = Arrays.copyOf(firstReadCountsByCycle, length);
                secondReadTotalsByCycle = Arrays.copyOf(secondReadTotalsByCycle, length);
                secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
                maxLengthSoFar = length;
            }
        }

        Histogram<Integer> getMeanQualityHistogram() {
            final String label = useOriginalQualities ? "MEAN_ORIGINAL_QUALITY" : "MEAN_QUALITY";
            final Histogram<Integer> meanQualities = new Histogram<>("CYCLE", label);

            int firstReadLength = 0;

            for (int cycle=0; cycle < firstReadTotalsByCycle.length; ++cycle) {
                if (firstReadTotalsByCycle[cycle] > 0) {
                    meanQualities.increment(cycle, firstReadTotalsByCycle[cycle] / firstReadCountsByCycle[cycle]);
                    firstReadLength = cycle;
                }
            }

            for (int i=0; i< secondReadTotalsByCycle.length; ++i) {
                if (secondReadCountsByCycle[i] > 0) {
                    final int cycle = firstReadLength + i;
                    meanQualities.increment(cycle, secondReadTotalsByCycle[i] / secondReadCountsByCycle[i]);
                }
            }

            return meanQualities;
        }

        boolean isEmpty() {
            return maxLengthSoFar == 0;
        }
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        if (PRODUCE_PLOT) {
            IOUtil.assertFileIsWritable(CHART_OUTPUT);
        }
        // If we're working with a single library, assign that library's name
        // as a suffix to the plot title
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // Skip unwanted records
        if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) return;
        if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) return;
        if (rec.isSecondaryOrSupplementary()) return;

        q.addRecord(rec);
        oq.addRecord(rec);
    }

    @Override
    protected void finish() {
        // Generate a "Histogram" of mean quality and write it to the file
        final MetricsFile<?,Integer> metrics = getMetricsFile();
        metrics.addHistogram(q.getMeanQualityHistogram());
        if (!oq.isEmpty()) metrics.addHistogram(oq.getMeanQualityHistogram());
        metrics.write(OUTPUT);

        if (q.isEmpty() && oq.isEmpty()) {
            logger.warn("No valid bases found in input file. No plot will be produced.");
        }
        else if (PRODUCE_PLOT){
            // Now run R to generate a chart
            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(new Resource(R_SCRIPT, MeanQualityByCycle.class));
            executor.addArgs(OUTPUT.getAbsolutePath(), CHART_OUTPUT.getAbsolutePath(), INPUT.getName(), plotSubtitle);
            executor.exec();
        }
    }
}

