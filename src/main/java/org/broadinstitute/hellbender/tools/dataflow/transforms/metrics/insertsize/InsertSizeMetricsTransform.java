package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.insertsize;

import com.google.cloud.dataflow.sdk.coders.BigEndianIntegerCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.Combine.CombineFn;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.Filter;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.cloud.genomics.dataflow.coders.GenericJsonCoder;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.HistogramCombinerDataflow;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.HistogramDataflow;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.MetricsFileDataflow;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

public class InsertSizeMetricsTransform extends PTransform<PCollection<GATKRead>, PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>>> {
    public static final long serialVersionUID = 1l;

    private final Arguments args;
    private final PCollectionView<List<Header>> metricHeaders;
    private final PCollectionView<SAMFileHeader> samHeader;

    /**
     * This selects the second in pair in order to avoid double counting inserts from pairs.
     * Ignore any non pair, duplicate or 0 length fragment because they don't make sense for count insert size.
     */
    private static final ReadFilter isSecondInMappedPair = r -> r.isPaired() &&
            !r.isUnmapped() &&
            !r.mateIsUnmapped() &&
            !r.isFirstOfPair() &&
            !(r.isSupplementaryAlignment() || r.isSecondaryAlignment()) &&
            !r.isDuplicate() &&
            r.getFragmentLength() != 0;

    public static class Arguments implements ArgumentCollectionDefinition, Serializable {
        public final static long serialVersionUID = 1l;

        @Argument(doc = "Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. " +
                "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
                "artifacts to make the mean and sd grossly misleading regarding the real distribution.")
        public double DEVIATIONS = 10;

        @Argument(shortName = "W", doc = "Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. " +
                "Also, when calculating mean and standard deviation, only bins <= HISTOGRAM_WIDTH will be included.", optional = true)
        public Integer HISTOGRAM_WIDTH = null;

        @Argument(shortName = "M", doc = "When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this " +
                "percentage of overall reads. (Range: 0 to 0.5).")
        public float MINIMUM_PCT = 0.05f;

        @Argument(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.  ")
        public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = EnumSet.of(MetricAccumulationLevel.ALL_READS);

        @Override
        public void validate() {
            if (MINIMUM_PCT < 0 || MINIMUM_PCT > 0.5) {
                throw new UserException.BadArgumentValue("MINIMUM_PCT", "It must be between 0 and 0.5 so all data categories don't get discarded.");
            }
        }
    }

    public InsertSizeMetricsTransform(Arguments args, PCollectionView<SAMFileHeader> samHeader, PCollectionView<List<Header>> metricHeaders) {
        args.validate();
        this.args = args;
        this.samHeader = samHeader;
        this.metricHeaders = metricHeaders;
    }

    @Override
    public PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> apply(PCollection<GATKRead> input) {
        PCollection<GATKRead> filtered = input.apply("Filter singletons and first of pair", Filter.by(isSecondInMappedPair));

        PCollection<KV<InsertSizeAggregationLevel, Integer>> kvPairs = filtered.apply(ParDo.named("Calculate metric and key")
                .withSideInputs(samHeader)
                .of(new DoFn<GATKRead, KV<InsertSizeAggregationLevel, Integer>>() {
                    private final static long serialVersionUID = 1l;

                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        final GATKRead read = c.element();
                        final SAMFileHeader samFileHeader = c.sideInput(samHeader);
                        Integer metric = computeMetric(read);
                        List<InsertSizeAggregationLevel> aggregationLevels = InsertSizeAggregationLevel.getKeysForAllAggregationLevels(read, samFileHeader, args.METRIC_ACCUMULATION_LEVEL);

                        aggregationLevels.stream().forEach(k -> c.output(KV.of(k, metric)));
                    }
                })).setCoder(KvCoder.of(GenericJsonCoder.of(InsertSizeAggregationLevel.class), BigEndianIntegerCoder.of()));

        CombineFn<Integer, HistogramDataflow<Integer>, HistogramDataflow<Integer>> combiner = new HistogramCombinerDataflow<>();
        PCollection<KV<InsertSizeAggregationLevel, HistogramDataflow<Integer>>> histograms = kvPairs.apply("Add reads to histograms",
                Combine.<InsertSizeAggregationLevel, Integer, HistogramDataflow<Integer>>perKey(combiner));

        PCollection<KV<InsertSizeAggregationLevel, KV<InsertSizeAggregationLevel, HistogramDataflow<Integer>>>> reKeyedHistograms = histograms.apply(ParDo.named("Re-key histograms").
                of(new DoFn<KV<InsertSizeAggregationLevel, HistogramDataflow<Integer>>, KV<InsertSizeAggregationLevel, KV<InsertSizeAggregationLevel, HistogramDataflow<Integer>>>>() {
                    public final static long serialVersionUID = 1l;

                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        KV<InsertSizeAggregationLevel, HistogramDataflow<Integer>> histo = c.element();
                        InsertSizeAggregationLevel oldKey = histo.getKey();
                        InsertSizeAggregationLevel newKey = new InsertSizeAggregationLevel(null, oldKey.getLibrary(), oldKey.getReadGroup(), oldKey.getSample());
                        c.output(KV.of(newKey, histo));
                    }
                }));

        PCollection<KV<InsertSizeAggregationLevel, MetricsFileDataflow<InsertSizeMetrics, Integer>>> metricsFiles = reKeyedHistograms.apply("Add histograms and metrics to MetricsFile",
                Combine.perKey(new CombineInsertSizeHistogramsIntoMetricsFile(args.DEVIATIONS, args.HISTOGRAM_WIDTH, args.MINIMUM_PCT)));

        PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> metricsFilesNoKeys = metricsFiles.apply(ParDo.named("Drop keys")
                .of(new DoFn<KV<?, MetricsFileDataflow<InsertSizeMetrics, Integer>>, MetricsFileDataflow<InsertSizeMetrics, Integer>>() {
                    public final static long serialVersionUID = 1l;

                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        c.output(c.element().getValue());
                    }
                }));

        PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> singleMetricsFile = metricsFilesNoKeys.<PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>>>apply("Combine aggregation levels",
                Combine.globally(new CombineInsertSizeMetricsFiles()));

        PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> singleMetricsFileWithHeaders = singleMetricsFile.apply(
                ParDo.named("add headers to MetricsFile")
                        .withSideInputs(this.metricHeaders)
                        .of(new AddHeadersToMetricsFile(metricHeaders)));

        return singleMetricsFileWithHeaders;
    }

    private Integer computeMetric(GATKRead read) {
        return Math.abs(read.getFragmentLength());
    }

    private static class AddHeadersToMetricsFile extends DoFn<MetricsFileDataflow<InsertSizeMetrics, Integer>, MetricsFileDataflow<InsertSizeMetrics, Integer>> {
        private final static long serialVersionUID = 1l;
        private PCollectionView<List<Header>> metricHeaders;

        private AddHeadersToMetricsFile(PCollectionView<List<Header>> metricHeaders) {
            this.metricHeaders = metricHeaders;
        }

        @Override
        public void processElement(ProcessContext c) throws Exception {
            final List<Header> headers = c.sideInput(metricHeaders);
            final MetricsFileDataflow<InsertSizeMetrics, Integer> metricFile = c.element();
            headers.stream().forEach(metricFile::addHeader);
            c.output(metricFile);
        }
    }
}
