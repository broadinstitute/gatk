package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
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
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.Serializable;
import java.util.Collections;
import java.util.List;

/**
 * Charts quality score distribution within a BAM file.
 *
 * This is the spark version of the tool.
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Program to chart " +
                "quality score distributions in a SAM/BAM file.",
        oneLineSummary = "QualityScoreDistribution on Spark",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@BetaFeature
public final class QualityScoreDistributionSpark extends GATKSparkTool {
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

    @Argument(shortName = "NC", fullName = "includeNoCalls", doc="If set to true, include quality for no-call bases in the distribution.")
    public boolean includeNoCalls = false;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @VisibleForTesting
    static final class Counts implements Serializable {
        private static final long serialVersionUID = 1L;
        public static final int MAX_BASE_QUALITY = 127;
        private final long[] qCounts;
        private final long[] oqCounts;
        private final boolean includeNoCalls; //NoCalls are 'N' (or 'n') bases in reads

        Counts(final boolean includeNoCalls){
           qCounts = new long[MAX_BASE_QUALITY + 1];
           oqCounts = new long[MAX_BASE_QUALITY + 1];
           this.includeNoCalls = includeNoCalls;
        }

        /**
         * Adds a read to this count object by increasing counts of the base qualities.
         */
        Counts addRead(final GATKRead read) {
            final byte[] bases = read.getBases();
            final byte[] quals = read.getBaseQualities();
            final byte[] oq    = ReadUtils.getOriginalBaseQualities(read);

            final int length = quals.length;

            for (int i=0; i<length; ++i) {
                if (includeNoCalls || !SequenceUtil.isNoCall(bases[i])) {
                    qCounts[quals[i]]++;
                    if (oq != null) {
                        oqCounts[oq[i]]++;
                    }
                }
            }
            return this;
        }

        /**
         * Adds the other count object into this one by adding the base quality counts.
         */
        Counts merge(final Counts counts2) {
            for (int i = 0; i <= MAX_BASE_QUALITY; i++) {
                this.qCounts[i] += counts2.qCounts[i];
            }
            for (int i = 0; i <= MAX_BASE_QUALITY; i++) {
                this.oqCounts[i] += counts2.oqCounts[i];
            }
            return this;
        }

        long[] getQualCounts(){
            return qCounts;
        }
        long[] getOrigQualCounts(){
            return oqCounts;
        }
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> reads = getReads();
        final MetricsReadFilter metricsFilter =
            new MetricsReadFilter(this.pfReadsOnly, this.alignedReadsOnly);
        final JavaRDD<GATKRead> filteredReads = reads.filter(read -> metricsFilter.test(read));
        final Counts result = filteredReads.aggregate(new Counts(includeNoCalls),
                (counts, read) -> counts.addRead(read),
                (counts1, counts2) -> counts1.merge(counts2));

        final MetricsFile<?, Byte> metrics = makeMetrics(result);
        saveResults(metrics, getHeaderForReads(), getReadSourceName());
    }

    //Convert the count object into a metrics object so save in a report
    private MetricsFile<?, Byte> makeMetrics(final Counts result) {
        // Built the Histograms out of the long[]s
        final Histogram<Byte> qHisto  = new Histogram<>("QUALITY", "COUNT_OF_Q");
        final Histogram<Byte> oqHisto = new Histogram<>("QUALITY", "COUNT_OF_OQ");

        for (int i=0; i< result.qCounts.length; ++i) {
            //Note: the checks are done to avoid work if there's nothing to add to the histogram
            if (result.qCounts[i]  > 0) qHisto.increment( (byte) i, (double) result.qCounts[i]);
            if (result.oqCounts[i] > 0) oqHisto.increment((byte) i, (double) result.oqCounts[i]);
        }

        final MetricsFile<?,Byte> metrics = getMetricsFile();
        metrics.addHistogram(qHisto);
        if (!oqHisto.isEmpty()) {
            metrics.addHistogram(oqHisto);
        }
        return metrics;
    }

    private void saveResults(final MetricsFile<?, Byte> metrics,  final SAMFileHeader readsHeader, final String inputFileName) {
        MetricsUtils.saveMetrics(metrics, out);

        if (metrics.getAllHistograms().isEmpty()) {
            logger.warn("No valid bases found in input file.");
        } else if(chartOutput != null){
            // If we're working with a single library, assign that library's name
            // as a suffix to the plot title
            String plotSubtitle = "";
            final List<SAMReadGroupRecord> readGroups = readsHeader.getReadGroups();
            if (readGroups.size() == 1) {
                plotSubtitle = readGroups.get(0).getLibrary();
                if (null == plotSubtitle) {
                    plotSubtitle = "";
                }
            }

            // Now run R to generate a chart
            //out is the metrics file
            //chartOutput is the pdf file with the graph
            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(getQualityScoreDistributionRScriptResource());
            executor.addArgs(out, chartOutput.getAbsolutePath(), inputFileName, plotSubtitle);
            executor.exec();
        }
    }

    /**
     * Return the R Script resource used by this class
     */
    @VisibleForTesting
    static Resource getQualityScoreDistributionRScriptResource() {
        final String R_SCRIPT = "qualityScoreDistribution.R";
        return new Resource(R_SCRIPT, QualityScoreDistributionSpark.class);
    }

}
