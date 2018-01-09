package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.StringUtil;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.filters.MetricsReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.metrics.MetricsUtils;
import org.broadinstitute.hellbender.metrics.analysis.BaseDistributionByCycleMetrics;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

/**
 * Collects base distribution per cycle in SAM/BAM/CRAM file(s). The tool leverages the Spark framework for
 * faster operation.
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CollectBaseDistributionByCycleSpark \
 *   -R gs://cloud-bucket/reference.fasta \
 *   -I gs://cloud-bucket/input.bam \
 *   -O gs://cloud-bucket/metrics_file.txt \
 *   -- \
 *   --spark-runner GCS \
 *   --cluster my-dataproc-cluster
 * </pre>
 *
 * See <a href ="https://software.broadinstitute.org/gatk/documentation/article?id=10060">Tutorial#10060</a>
 * for an example of how to set up and run a Spark tool on a cloud Spark cluster.
 */

@CommandLineProgramProperties(
        summary = "Collects base distribution per cycle in SAM/BAM/CRAM file(s). The tool " +
                "leverages the Spark framework for faster operation.",
        oneLineSummary = "Collects base distribution per cycle in SAM/BAM/CRAM file(s).",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class CollectBaseDistributionByCycleSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(
            doc = "Output metrics file.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true
    )
    public String out;

    @Argument(
            shortName="C",
            fullName = "chart",
            doc="Output charts file (pdf).",
            optional=true
    )
    public File chartOutput;

    @Argument(
            shortName="A",
            fullName = "aligned-reads-only",
            doc="If set to true, calculates the base distribution over aligned reads only."
    )
    public boolean alignedReadsOnly = false;

    @Argument(
            shortName="F",
            fullName = "pf-reads-only",
            doc="If set to true, calculates the base distribution over passing filter (PF) reads only."
    )
    public boolean pfReadsOnly = false;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Arrays.asList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    private static final class HistogramGenerator implements Serializable{
        private static final long serialVersionUID = 1L;

        private static final int NBASES = 5; //ACGT + unknown
        private int maxLengthSoFar = 0;
        private final long[][] firstReadTotalsByCycle = new long[NBASES][maxLengthSoFar];
        private long[] firstReadCountsByCycle = new long[maxLengthSoFar];
        private final long[][] secondReadTotalsByCycle = new long[NBASES][maxLengthSoFar];
        private long[] secondReadCountsByCycle = new long[maxLengthSoFar];
        private boolean seenAnySecondOfPair = false;

        /**
         * Merges the argument into this histogram generator. Returns the modified 'this' object.
         */
        public HistogramGenerator merge(final HistogramGenerator hg2) {
            Utils.nonNull(hg2);
            ensureArraysBigEnough(hg2.maxLengthSoFar);
            for (int i = 0; i < NBASES; i++) {
                for (int j = 0; j < hg2.firstReadTotalsByCycle[i].length; j++) {
                    this.firstReadTotalsByCycle[i][j] += hg2.firstReadTotalsByCycle[i][j];
                }
                for (int j = 0; j < hg2.secondReadTotalsByCycle[i].length; j++) {
                    this.secondReadTotalsByCycle[i][j] += hg2.secondReadTotalsByCycle[i][j];
                }
            }
            for (int i = 0; i < hg2.firstReadCountsByCycle.length; i++) {
                this.firstReadCountsByCycle[i] += hg2.firstReadCountsByCycle[i];
            }
            for (int i = 0; i < hg2.secondReadCountsByCycle.length; i++) {
                this.secondReadCountsByCycle[i] += hg2.secondReadCountsByCycle[i];
            }
            this.seenAnySecondOfPair = this.seenAnySecondOfPair || hg2.seenAnySecondOfPair;
            return this;
        }

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

        /**
         * Adds one read to this histogram.
         */
        HistogramGenerator addRead(final GATKRead rec) {
            final byte[] bases = rec.getBases();
            if (bases == null) {
                return this;
            }
            final int length = bases.length;
            final boolean rc = rec.isReverseStrand();
            ensureArraysBigEnough(length + 1);
            if (rec.isPaired() && rec.isSecondOfPair()) {
                seenAnySecondOfPair = true;
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
            return this;
        }

        private void ensureArraysBigEnough(final int length) {
            if (length > maxLengthSoFar) {
                for (int i = 0; i < NBASES; i++) {
                    firstReadTotalsByCycle[i] = Arrays.copyOf(firstReadTotalsByCycle[i], length);
                    secondReadTotalsByCycle[i] = Arrays.copyOf(secondReadTotalsByCycle[i], length);
                }
                firstReadCountsByCycle = Arrays.copyOf(firstReadCountsByCycle, length);
                secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
                maxLengthSoFar = length;
            }
        }

        public void addToMetricsFile(final MetricsFile<BaseDistributionByCycleMetrics, ?> metrics) {
            int maxFirstReadCycles = 0;
            for (int i = 0; i < maxLengthSoFar; i++) {
                if (firstReadCountsByCycle[i] != 0) {
                    final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                    metric.READ_END = 1;
                    metric.CYCLE = i;
                    metric.PCT_A = (100.0 * firstReadTotalsByCycle[0][i] / firstReadCountsByCycle[i]);
                    metric.PCT_C = (100.0 * firstReadTotalsByCycle[1][i] / firstReadCountsByCycle[i]);
                    metric.PCT_G = (100.0 * firstReadTotalsByCycle[2][i] / firstReadCountsByCycle[i]);
                    metric.PCT_T = (100.0 * firstReadTotalsByCycle[3][i] / firstReadCountsByCycle[i]);
                    metric.PCT_N = (100.0 * firstReadTotalsByCycle[4][i] / firstReadCountsByCycle[i]);
                    metrics.addMetric(metric);
                    maxFirstReadCycles = i;
                }
            }
            if (seenAnySecondOfPair) {
                for (int i = 0; i < maxLengthSoFar; i++) {
                    if (secondReadCountsByCycle[i] != 0) {
                        final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                        metric.READ_END = 2;
                        metric.CYCLE = (i + maxFirstReadCycles);
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

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        final MetricsFile<BaseDistributionByCycleMetrics, Integer> metricsFile = calculateBaseDistributionByCycle(reads);
        saveResults(metricsFile, getHeaderForReads(), getReadSourceName());
    }

    /**
     * Computes the MeanQualityByCycle. Creates a metrics file with relevant histograms.
     */
    public MetricsFile<BaseDistributionByCycleMetrics, Integer> calculateBaseDistributionByCycle(final JavaRDD<GATKRead> reads){
        final MetricsReadFilter metricsFilter =
            new MetricsReadFilter(this.pfReadsOnly, this.alignedReadsOnly);
        final JavaRDD<GATKRead> filteredReads = reads.filter(read -> metricsFilter.test(read));
        final HistogramGenerator hist = filteredReads.aggregate(new HistogramGenerator(),
                (hgp, read) -> hgp.addRead(read),
                (hgp1, hgp2) -> hgp1.merge(hgp2));

        final MetricsFile<BaseDistributionByCycleMetrics, Integer> metricsFile = getMetricsFile();
        hist.addToMetricsFile(metricsFile);
        return metricsFile;
    }

    protected void saveResults(final MetricsFile<?, Integer> metrics, final SAMFileHeader readsHeader, final String inputFileName) {
        MetricsUtils.saveMetrics(metrics, out);

        if (metrics.getAllHistograms().isEmpty()) {
            logger.warn("No valid bases found in input file.");
        } else if (chartOutput != null) {
            // Now run R to generate a chart

            // If we're working with a single library, assign that library's name
            // as a suffix to the plot title
            final List<SAMReadGroupRecord> readGroups = readsHeader.getReadGroups();

            /*
             * A subtitle for the plot, usually corresponding to a library.
             */
            String plotSubtitle = "";
            if (readGroups.size() == 1) {
                plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
            }
            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(getBaseDistributionByCycleRScriptResource());
            executor.addArgs(out, chartOutput.getAbsolutePath(), inputFileName, plotSubtitle);
            executor.exec();
        }
    }

    /**
     * Return the R Script resource used by this class
     */
    @VisibleForTesting
    static Resource getBaseDistributionByCycleRScriptResource() {
        final String R_SCRIPT = "baseDistributionByCycle.R";
        return new Resource(R_SCRIPT, CollectBaseDistributionByCycleSpark.class);
    }

}
