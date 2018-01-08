package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
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
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Program to generate a data table and chart of mean quality by cycle from a
 * BAM file.  Works best on a single lane/run of data, but can be applied to
 * merged BAMs - the output may just be a little confusing.
 *
 * This is the Spark implementation of this tool.
 */
@DocumentedFeature
@CommandLineProgramProperties(
summary = "Program to generate a data table and pdf chart of " +
        "mean base quality by cycle from a SAM/BAM file.  Works best on a single lane/run of data, but can be applied to " +
        "merged BAMs. Uses R to generate chart output.",
        oneLineSummary = "MeanQualityByCycle on Spark",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@BetaFeature
public final class MeanQualityByCycleSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1l;

    @Argument(doc = "uri for the output file: a local file path",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public String out;

    @Argument(shortName="C", fullName = "chart", doc="A file (with .pdf extension) to write the chart to.", optional=true)
    public File chartOutput;

    @Argument(shortName="A", fullName = "alignedReadsOnly", doc="If set to true calculate mean quality over aligned reads only.")
    public boolean alignedReadsOnly = false;

    @Argument(shortName="F", fullName = "pfReadsOnly", doc="If set to true calculate mean quality over PF reads only.")
    public boolean pfReadsOnly = false;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @VisibleForTesting
    static final class HistogramGenerator implements Serializable {
        private static final long serialVersionUID = 1L;
        final boolean useOriginalQualities;
        int maxLengthSoFar = 0;
        double[] firstReadTotalsByCycle  = new double[maxLengthSoFar];
        long[]   firstReadCountsByCycle  = new long[maxLengthSoFar];
        double[] secondReadTotalsByCycle = new double[maxLengthSoFar];
        long[]   secondReadCountsByCycle = new long[maxLengthSoFar];

        HistogramGenerator(final boolean useOriginalQualities) {
            this.useOriginalQualities = useOriginalQualities;
        }

        HistogramGenerator addRead(final GATKRead read) {
            final byte[] quals = (useOriginalQualities ? ReadUtils.getOriginalBaseQualities(read) : read.getBaseQualities());
            if (quals == null) {
                return this;
            }

            final int length = quals.length;
            final boolean isReverseStrand = read.isReverseStrand();
            ensureArraysBigEnough(length+1);

            for (int i=0; i<length; ++i) {
                final int cycle = isReverseStrand ? length-i : i+1;

                if (read.isPaired() && read.isSecondOfPair()) {
                    secondReadTotalsByCycle[cycle] += quals[i];
                    secondReadCountsByCycle[cycle] += 1;
                }
                else {
                    firstReadTotalsByCycle[cycle] += quals[i];
                    firstReadCountsByCycle[cycle] += 1;
                }
            }
            return this;
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


        /**
         * Merges the argument into this histogram generator. Returns the modified 'this' object.
         * Note: you can only merge HistogramGenerator is they have the same 'useOriginalQualities' value.
         */
        public HistogramGenerator merge(final HistogramGenerator hg2) {
            Utils.nonNull(hg2);
            Utils.validateArg(this.useOriginalQualities == hg2.useOriginalQualities,
                    () -> "unequal useOriginalQualities. This has " + this.useOriginalQualities);

            ensureArraysBigEnough(hg2.maxLengthSoFar);
            for (int i = 0; i < hg2.firstReadTotalsByCycle.length; i++) {
               this.firstReadTotalsByCycle[i] += hg2.firstReadTotalsByCycle[i];
            }
            for (int i = 0; i < hg2.secondReadTotalsByCycle.length; i++) {
                this.secondReadTotalsByCycle[i] += hg2.secondReadTotalsByCycle[i];
            }
            for (int i = 0; i < hg2.firstReadCountsByCycle.length; i++) {
                this.firstReadCountsByCycle[i] += hg2.firstReadCountsByCycle[i];
            }
            for (int i = 0; i < hg2.secondReadCountsByCycle.length; i++) {
                this.secondReadCountsByCycle[i] += hg2.secondReadCountsByCycle[i];
            }
            return this;
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

    /**
     * This pair makes it easy to process the set of reads once and maintain both histograms.
     */
    private static final class HistogramGeneratorPair implements Serializable {
        private static final long serialVersionUID = 1;
        HistogramGenerator useQuals = new HistogramGenerator(false);
        HistogramGenerator useOrigQuals = new HistogramGenerator(true);
        HistogramGeneratorPair addRead(final GATKRead read){
            useQuals.addRead(read);
            useOrigQuals.addRead(read);
            return this;
        }
        HistogramGeneratorPair merge(final HistogramGeneratorPair other){
            useQuals.merge(other.useQuals);
            useOrigQuals.merge(other.useOrigQuals);
            return this;
        }
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        final MetricsFile<?, Integer> metricsFile = calculateMeanQualityByCycle(reads);
        saveResults(metricsFile, getHeaderForReads(), getReadSourceName());
    }

    /**
     * Computes the MeanQualityByCycle. Creates a metrics file with relevant histograms.
     */
    public MetricsFile<?, Integer> calculateMeanQualityByCycle(final JavaRDD<GATKRead> reads){
        final MetricsReadFilter metricsFilter =
            new MetricsReadFilter(this.pfReadsOnly, this.alignedReadsOnly);
        final JavaRDD<GATKRead> filteredReads = reads.filter(read -> metricsFilter.test(read));
        final HistogramGeneratorPair aggregate = filteredReads.aggregate(new HistogramGeneratorPair(),
                (hgp, read) -> hgp.addRead(read),
                (hgp1, hgp2) -> hgp1.merge(hgp2));
        return finish(aggregate.useQuals, aggregate.useOrigQuals);
    }

    private MetricsFile<?,Integer> finish(final HistogramGenerator q, final HistogramGenerator oq) {
        // Generate a "Histogram" of mean quality and write it to the file
        final MetricsFile<?,Integer> metrics = getMetricsFile();
        metrics.addHistogram(q.getMeanQualityHistogram());
        if (!oq.isEmpty()) {
            metrics.addHistogram(oq.getMeanQualityHistogram());
        }
        return metrics;
    }

    private void saveResults(final MetricsFile<?, Integer> metrics, final SAMFileHeader readsHeader, final String inputFileName){
        MetricsUtils.saveMetrics(metrics, out);

        if (metrics.getAllHistograms().isEmpty()) {
            logger.warn("No valid bases found in input file.");
        } else if (chartOutput != null){
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
            executor.addScript(getMeanQualityByCycleRScriptResource());
            executor.addArgs(out, chartOutput.getAbsolutePath(), inputFileName, plotSubtitle);
            executor.exec();
        }
    }

    /**
     * Return the R Script resource used by this class
     */
    @VisibleForTesting
    static Resource getMeanQualityByCycleRScriptResource() {
        final String R_SCRIPT = "meanQualityByCycle.R";
        return new Resource(R_SCRIPT, MeanQualityByCycleSpark.class);
    }

}
