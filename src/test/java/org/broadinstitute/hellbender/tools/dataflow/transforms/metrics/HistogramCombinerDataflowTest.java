package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics;

import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.util.SerializableUtils;
import htsjdk.samtools.util.TestUtil;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class HistogramCombinerDataflowTest {

    @Test
    public void testHistogramCombiner(){
        List<Integer> records = IntStream.rangeClosed(1, 1000).boxed().collect(Collectors.toList());
        Combine.CombineFn<Integer, HistogramDataflow<Integer>, HistogramDataflow<Integer>> combiner = new HistogramCombinerDataflow<>();

        HistogramDataflow<Integer> result = combiner.apply(records);
        Assert.assertEquals(result.getCount(), 1000.0);
        Assert.assertEquals(result.getMax(),1000.0);
        Assert.assertEquals(result.getMin(),1.0);
    }


    @Test
    public void dataflowSerializeMetricsFileTest(){
        MetricsFileDataflow<InsertSizeMetrics,Integer> metrics = new MetricsFileDataflow<>();
        metrics.addHistogram(new HistogramDataflow<>());

        @SuppressWarnings("unchecked")
        MetricsFileDataflow<InsertSizeMetrics, Integer> newMetrics =
                SerializableUtils.ensureSerializableByCoder(SerializableCoder.of(MetricsFileDataflow.class), metrics, "error");
        Assert.assertEquals(newMetrics.getAllHistograms(),metrics.getAllHistograms());
    }

    @Test
    public void javaSerializeMetricsFileTest() throws IOException, ClassNotFoundException {
        final MetricsFileDataflow<InsertSizeMetrics,Integer> metrics = new MetricsFileDataflow<>();
        metrics.addHistogram(new HistogramDataflow<>());
        final MetricsFileDataflow<InsertSizeMetrics,Integer> deserializedMetrics = TestUtil.serializeAndDeserialize(metrics);

        Assert.assertEquals(deserializedMetrics.getAllHistograms(), metrics.getAllHistograms());
    }
}