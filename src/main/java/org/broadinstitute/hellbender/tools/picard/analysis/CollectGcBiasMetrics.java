package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.text.NumberFormat;
import java.util.*;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Tool to collect information about GC bias in the reads in a given SAM/BAM file. Computes" +
                " the number of windows (of size specified by WINDOW_SIZE) in the genome at each GC%" +
                " and counts the number of read starts in each GC bin.  What is output and plotted is" +
                " the \"normalized coverage\" in each bin - i.e. the number of reads per window normalized" +
                " to the average number of reads per window across the whole genome..\n",
        oneLineSummary = "Produces metrics about GC bias in the reads in the provided SAM/BAM file",
        programGroup = QCProgramGroup.class
)
public final class CollectGcBiasMetrics extends SinglePassSamProgram {
    /** The location of the R script to do the plotting. */
    private static final String R_SCRIPT = "gcBias.R";

    // Usage and parameters

    @Argument(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Argument(shortName = "S", doc = "The text file to write summary metrics to.", optional = true)
    public File SUMMARY_OUTPUT;

    @Argument(doc = "The size of windows on the genome that are used to bin reads.")
    public int WINDOW_SIZE = 100;

    @Argument(doc = "For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    @Argument(shortName = "BS", doc = "Whether the SAM/BAM file consists of bisulfite sequenced reads.  ")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Argument(doc = "Should an output plot be created")
    public boolean PRODUCE_PLOT = false;

    // Used to keep track of the total clusters as this is kinda important for bias
    private int totalClusters = 0;
    private int totalAlignedReads = 0;
    // Histograms to track the number of windows at each GC, and the number of read starts
    // at windows of each GC. Need 101 to get from 0-100.
    private static final int WINDOWS = 101;
    private final int[] windowsByGc = new int[WINDOWS];
    private final int[] readsByGc = new int[WINDOWS];
    private final long[] basesByGc = new long[WINDOWS];
    private final long[] errorsByGc = new long[WINDOWS];
    private static int lastContig = SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
    private byte[] gc;
    private byte[] refBases;
    private String saveHeader;

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        if (SUMMARY_OUTPUT != null) IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);
        saveHeader = header.getReadGroups().get(0).getLibrary();
    }

    ////////////////////////////////////////////////////////////////////////////
    // Loop over the reference and the reads and calculate the basic metrics
    ////////////////////////////////////////////////////////////////////////////
    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        //if read is unaligned then ref is passed in as null
        if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++this.totalClusters;
        if (ref!=null) {
            //only do the recalculation of gc if current ref is different from last ref
            if (ref.getContigIndex() != lastContig) {
                refBases = ref.getBases();
                StringUtil.toUpperCase(refBases);
                final int refLength = refBases.length;
                final int lastWindowStart = refLength - WINDOW_SIZE;
                gc = calculateAllGcs(refBases, windowsByGc, lastWindowStart);
                lastContig = ref.getContigIndex();
            }
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) ++this.totalClusters;
            if (!rec.getReadUnmappedFlag()) {
                final int pos = rec.getReadNegativeStrandFlag() ? rec.getAlignmentEnd() - WINDOW_SIZE : rec.getAlignmentStart();
                ++this.totalAlignedReads;
                if (pos > 0) {
                    final int windowGc = gc[pos];
                    if (windowGc >= 0) {
                        ++readsByGc[windowGc];
                        basesByGc[windowGc] += rec.getReadLength();
                        errorsByGc[windowGc] +=
                                SequenceUtil.countMismatches(rec, refBases, IS_BISULFITE_SEQUENCED) +
                                        SequenceUtil.countInsertedBases(rec) + SequenceUtil.countDeletedBases(rec);
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    // Synthesize the normalized coverage metrics and write it all out to a file
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void finish () {
        final MetricsFile<GcBiasDetailMetrics, ?> metricsFile = getMetricsFile();
        final double totalWindows = sum(windowsByGc);
        final double totalReads = sum(readsByGc);
        final double meanReadsPerWindow = totalReads / totalWindows;
        final double minimumWindowsToCountInSummary = totalWindows * this.MINIMUM_GENOME_FRACTION;

        for (int i = 0; i < windowsByGc.length; ++i) {
            if (windowsByGc[i] == 0) continue;

            final GcBiasDetailMetrics m = new GcBiasDetailMetrics();
            m.GC = i;
            m.WINDOWS = windowsByGc[i];
            m.READ_STARTS = readsByGc[i];
            if (errorsByGc[i] > 0) m.MEAN_BASE_QUALITY = QualityUtil.getPhredScoreFromObsAndErrors(basesByGc[i], errorsByGc[i]);
            m.NORMALIZED_COVERAGE = (m.READ_STARTS / (double) m.WINDOWS) / meanReadsPerWindow;
            m.ERROR_BAR_WIDTH = (Math.sqrt(m.READ_STARTS) / (double) m.WINDOWS) / meanReadsPerWindow;

            metricsFile.addMetric(m);
        }

        metricsFile.write(OUTPUT);

        // Synthesize the high level metrics
        if (SUMMARY_OUTPUT != null) {
            final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
            final GcBiasSummaryMetrics summary = new GcBiasSummaryMetrics();
            summary.WINDOW_SIZE = this.WINDOW_SIZE;
            summary.TOTAL_CLUSTERS = this.totalClusters;
            summary.ALIGNED_READS = this.totalAlignedReads;
            calculateDropoutMetrics(metricsFile.getMetrics(), summary);

            summaryMetricsFile.addMetric(summary);
            summaryMetricsFile.write(SUMMARY_OUTPUT);
        }

        // Plot the results
        final NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        final String subtitle = "Total clusters: " + fmt.format(this.totalClusters) +
                ", Aligned reads: " + fmt.format(this.totalAlignedReads);
        String title = INPUT.getName().replace(".duplicates_marked", "").replace(".aligned.bam", "");
        title += "." + saveHeader;
        if (PRODUCE_PLOT){
            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(new Resource(R_SCRIPT, CollectGcBiasMetrics.class));
            executor.addArgs(OUTPUT.getAbsolutePath(), CHART_OUTPUT.getAbsolutePath(), title, subtitle, String.valueOf(WINDOW_SIZE));
            executor.exec();
        }
    }

    /** Sums the values in an int[]. */
    private double sum(final int[] values) {
        final int length = values.length;
        double total = 0;
        for (int i = 0; i < length; ++i) {
            total += values[i];
        }

        return total;
    }

    /** Calculates the Illumina style AT and GC dropout numbers. */
    private void calculateDropoutMetrics(final Collection<GcBiasDetailMetrics> details,
                                         final GcBiasSummaryMetrics summary) {
        // First calculate the totals
        double totalReads = 0;
        double totalWindows = 0;

        for (final GcBiasDetailMetrics detail : details) {
            totalReads += detail.READ_STARTS;
            totalWindows += detail.WINDOWS;
        }

        double atDropout = 0;
        double gcDropout = 0;

        for (final GcBiasDetailMetrics detail : details) {
            final double relativeReads = detail.READ_STARTS / totalReads;
            final double relativeWindows = detail.WINDOWS / totalWindows;
            final double dropout = (relativeWindows - relativeReads) * 100;

            if (dropout > 0) {
                if (detail.GC <= 50) atDropout += dropout;
                if (detail.GC >= 50) gcDropout += dropout;
            }
        }

        summary.AT_DROPOUT = atDropout;
        summary.GC_DROPOUT = gcDropout;
    }

    /** Calculcate all the GC values for all windows. */
    private byte[] calculateAllGcs(final byte[] refBases, final int[] windowsByGc, final int lastWindowStart) {
        final int refLength = refBases.length;
        final byte[] gc = new byte[refLength + 1];
        final CalculateGcState state = new CalculateGcState();
        for (int i = 1; i < lastWindowStart; ++i) {
            final int windowEnd = i + WINDOW_SIZE;
            final int windowGc = calculateGc(refBases, i, windowEnd, state);
            gc[i] = (byte) windowGc;
            if (windowGc != -1) windowsByGc[windowGc]++;
        }
        return gc;
    }

    /**
     * Calculates GC as a number from 0 to 100 in the specified window. If the window includes
     * more than five no-calls then -1 is returned.
     */
    private int calculateGc(final byte[] bases, final int startIndex, final int endIndex, final CalculateGcState state) {
        if (state.init) {
            state.init = false;
            state.gcCount = 0;
            state.nCount = 0;
            for (int i = startIndex; i < endIndex; ++i) {
                final byte base = bases[i];
                if (base == 'G' || base == 'C') ++state.gcCount;
                else if (base == 'N') ++state.nCount;
            }
        } else {
            final byte newBase = bases[endIndex - 1];
            if (newBase == 'G' || newBase == 'C') ++state.gcCount;
            else if (newBase == 'N') ++state.nCount;

            if (state.priorBase == 'G' || state.priorBase == 'C') --state.gcCount;
            else if (state.priorBase == 'N') --state.nCount;
        }
        state.priorBase = bases[startIndex];
        if (state.nCount > 4) return -1;
        else return (state.gcCount * 100) / (endIndex - startIndex);
    }

    /** Keeps track of current GC calculation state. */
    class CalculateGcState {
        boolean init = true;
        int nCount;
        int gcCount;
        byte priorBase;
    }
}
