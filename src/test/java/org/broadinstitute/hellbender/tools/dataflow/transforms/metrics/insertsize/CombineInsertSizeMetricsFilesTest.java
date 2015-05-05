package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.insertsize;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Sets;
import htsjdk.samtools.metrics.StringHeader;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.HistogramDataflow;
import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.MetricsFileDataflow;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by louisb on 8/12/15.
 */
public class CombineInsertSizeMetricsFilesTest {

    @Test
    public void testCombineMetricsFiles(){
        Combine.CombineFn<MetricsFileDataflow<InsertSizeMetrics,Integer>,
                MetricsFileDataflow<InsertSizeMetrics,Integer>,
                MetricsFileDataflow<InsertSizeMetrics,Integer>> combiner =
                new CombineInsertSizeMetricsFiles();

        MetricsFileDataflow<InsertSizeMetrics, Integer> mf1 = createMetricFileWith1Metric("header1");
        MetricsFileDataflow<InsertSizeMetrics, Integer> mf2 = createMetricFileWith1Metric("header2");
        MetricsFileDataflow<InsertSizeMetrics, Integer> mf3 = createMetricFileWith1Metric("header3");

        @SuppressWarnings("unchecked")
        MetricsFileDataflow<InsertSizeMetrics, Integer> combined = combiner.apply(ImmutableList.of(mf1, mf2, mf3));

        Assert.assertEquals(combined.getMetrics().size(), 3);
        Assert.assertEquals(combined.getAllHistograms().size(), 0);
        Assert.assertEquals(combined.getHeaders().size(), 3);
    }

    @Test
    public void testCombineMetricsFilePTransform(){
        final Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        MetricsFileDataflow<InsertSizeMetrics, Integer> mf1 = createMetricFileWith1Metric("header1");
        mf1.addHistogram(createDummyHistogram());

        MetricsFileDataflow<InsertSizeMetrics, Integer> mf2 = createMetricFileWith1Metric("header2");
        mf2.addHistogram(createDummyHistogram());

        MetricsFileDataflow<InsertSizeMetrics, Integer> mf3 = createMetricFileWith1Metric("header3");

        PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> files = p.apply(Create.of(ImmutableList.of(mf1, mf2, mf3)));
        PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> combined = files.apply(Combine.globally(new CombineInsertSizeMetricsFiles()));
        DirectPipelineRunner.EvaluationResults results =  (DirectPipelineRunner.EvaluationResults)p.run();

        MetricsFileDataflow<InsertSizeMetrics, Integer> metricsFile = results.getPCollection(combined).get(0);

        Assert.assertEquals(Sets.newHashSet(metricsFile.getHeaders()), Sets.newHashSet(new StringHeader("header1"), new StringHeader("header2"), new StringHeader("header3")));
        Assert.assertEquals(metricsFile.getAllHistograms().size(), 2);
    }

    private static MetricsFileDataflow<InsertSizeMetrics, Integer> createMetricFileWith1Metric(String header) {
        MetricsFileDataflow<InsertSizeMetrics,Integer> mf1 = new MetricsFileDataflow<>();
        mf1.addMetric(new InsertSizeMetrics());
        mf1.addHeader(new StringHeader(header));
        return mf1;
    }

    private HistogramDataflow<Integer> createDummyHistogram() {
        HistogramDataflow<Integer> h1= new HistogramDataflow<>();
        h1.addInput(10);
        h1.addInput(20);
        h1.addInput(10);
        return h1;
    }
}