package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
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

@CommandLineProgramProperties(
        summary = "Program to chart the nucleotide distribution per cycle in a SAM/BAM/CRAM file",
        oneLineSummary = "Produces metrics about nucleotide distribution per cycle in a SAM/BAM/CRAM file",
        programGroup = QCProgramGroup.class
)
public final class CollectBaseDistributionByCycle extends SinglePassSamProgram {
    public static final String R_SCRIPT = "baseDistributionByCycle.R";

    @Argument(shortName = "CHART", doc = "A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Argument(doc = "If set to true, calculate the base distribution over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(doc = "If set to true calculate the base distribution over PF reads only.")
    public boolean PF_READS_ONLY = false;

    @Argument(doc = "Should an output plot be created")
    public boolean PRODUCE_PLOT = false;

    private HistogramGenerator hist;
    private String plotSubtitle = "";
    private final Logger log = LogManager.getLogger(CollectBaseDistributionByCycle.class);

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
        }
        hist = new HistogramGenerator();
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        if ((PF_READS_ONLY) && (rec.getReadFailsVendorQualityCheckFlag())) {
            return;
        }
        if ((ALIGNED_READS_ONLY) && (rec.getReadUnmappedFlag())) {
            return;
        }
        if (rec.isSecondaryOrSupplementary()) {
            return;
        }
        hist.addRecord(rec);
    }

    @Override
    protected void finish() {
        final MetricsFile<BaseDistributionByCycleMetrics, ?> metrics = getMetricsFile();
        hist.addToMetricsFile(metrics);
        metrics.write(OUTPUT);

        if (hist.isEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        } else if(PRODUCE_PLOT){
            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(new Resource(R_SCRIPT, CollectBaseDistributionByCycle.class));
            executor.addArgs(OUTPUT.getAbsolutePath(), CHART_OUTPUT.getAbsolutePath(), INPUT.getName(), plotSubtitle);
            executor.exec();
        }
    }

    private class HistogramGenerator {
        private int maxLengthSoFar = 0;
        private final long[][] firstReadTotalsByCycle = new long[5][maxLengthSoFar];
        private long[] firstReadCountsByCycle = new long[maxLengthSoFar];
        private final long[][] secondReadTotalsByCycle = new long[5][maxLengthSoFar];
        private long[] secondReadCountsByCycle = new long[maxLengthSoFar];
        private boolean seenSecondEnd = false;

        private int baseToInt(final byte base) {
            switch (base) {
                case 'A':
                case 'a':
                    return 0;
                case 'C':
                case 'c':
                    return 1;
                case 'G':
                case 'g':
                    return 2;
                case 'T':
                case 't':
                    return 3;
            }
            return 4;
        }

        void addRecord(final SAMRecord rec) {
            final byte[] bases = rec.getReadBases();
            if (bases == null) {
                return;
            }
            final int length = bases.length;
            final boolean rc = rec.getReadNegativeStrandFlag();
            ensureArraysBigEnough(length + 1);
            if ((rec.getReadPairedFlag()) && (rec.getSecondOfPairFlag())) {
                seenSecondEnd = true;
                for (int i = 0; i < length; i++) {
                    final int cycle = rc ? length - i : i + 1;
                    secondReadTotalsByCycle[baseToInt(bases[i])][cycle] += 1;
                    secondReadCountsByCycle[cycle] += 1;
                }
            } else {
                for (int i = 0; i < length; i++) {
                    final int cycle = rc ? length - i : i + 1;
                    firstReadTotalsByCycle[baseToInt(bases[i])][cycle] += 1;
                    firstReadCountsByCycle[cycle] += 1;
                }
            }
        }

        private void ensureArraysBigEnough(final int length) {
            if (length > maxLengthSoFar) {
                for (int i = 0; i < 5; i++) {
                    firstReadTotalsByCycle[i] = Arrays.copyOf(firstReadTotalsByCycle[i], length);
                    secondReadTotalsByCycle[i] = Arrays.copyOf(secondReadTotalsByCycle[i], length);
                }
                firstReadCountsByCycle = Arrays.copyOf(firstReadCountsByCycle, length);
                secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
                maxLengthSoFar = length;
            }
        }

        boolean isEmpty() {
            return maxLengthSoFar == 0;
        }

        public void addToMetricsFile(final MetricsFile<BaseDistributionByCycleMetrics, ?> metrics) {
            int firstReadLength = 0;
            for (int i = 0; i < maxLengthSoFar; i++) {
                if (0 != firstReadCountsByCycle[i]) {
                    final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                    metric.READ_END = 1;
                    metric.CYCLE = i;
                    metric.PCT_A = (100.0 * firstReadTotalsByCycle[0][i] / firstReadCountsByCycle[i]);
                    metric.PCT_C = (100.0 * firstReadTotalsByCycle[1][i] / firstReadCountsByCycle[i]);
                    metric.PCT_G = (100.0 * firstReadTotalsByCycle[2][i] / firstReadCountsByCycle[i]);
                    metric.PCT_T = (100.0 * firstReadTotalsByCycle[3][i] / firstReadCountsByCycle[i]);
                    metric.PCT_N = (100.0 * firstReadTotalsByCycle[4][i] / firstReadCountsByCycle[i]);
                    metrics.addMetric(metric);
                    firstReadLength = i;
                }
            }
            if (seenSecondEnd) {
                for (int i = 0; i < maxLengthSoFar; i++) {
                    if (0 != secondReadCountsByCycle[i]) {
                        final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                        metric.READ_END = 2;
                        metric.CYCLE = (i + firstReadLength);
                        metric.PCT_A = (100.0 * secondReadTotalsByCycle[0][i] / secondReadCountsByCycle[i]);
                        metric.PCT_C = (100.0 * secondReadTotalsByCycle[1][i] / secondReadCountsByCycle[i]);
                        metric.PCT_G = (100.0 * secondReadTotalsByCycle[2][i] / secondReadCountsByCycle[i]);
                        metric.PCT_T = (100.0 * secondReadTotalsByCycle[3][i] / secondReadCountsByCycle[i]);
                        metric.PCT_N = (100.0 * secondReadTotalsByCycle[4][i] / secondReadCountsByCycle[i]);
                        metrics.addMetric(metric);
                    }
                }
            }
        }
    }
}
