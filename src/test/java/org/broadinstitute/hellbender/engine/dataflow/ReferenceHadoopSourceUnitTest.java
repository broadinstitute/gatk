package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.VoidCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceHadoopSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.io.Serializable;
import java.util.Map;

public class ReferenceHadoopSourceUnitTest extends BaseTest implements Serializable {
    private static final long serialVersionUID = 1L;

    @Test
    public void testReferenceSourceQuery() throws IOException {
        // test ReferenceHadoopSource by reading a reference from a local file (not HDFS)
        final String referencePath = hg19MiniReference;

        final ReferenceHadoopSource refSource = new ReferenceHadoopSource(referencePath);
        final ReferenceBases bases = refSource.getReferenceBases(PipelineOptionsFactory.create(),
                new SimpleInterval("2", 10001, 10010));

        Assert.assertNotNull(bases);
        Assert.assertEquals(bases.getBases().length, 10, "Wrong number of bases returned");
        Assert.assertEquals(new String(bases.getBases()), "CGTATCCCAC", "Wrong bases returned");
    }

    @Test
    public void testBroadcastReference() throws IOException {
        final String referencePath = hg19MiniReference;
        final Pipeline p = GATKTestPipeline.create();

        final ReferenceFileSource refSource = new ReferenceFileSource(referencePath);
        PCollection<KV<String, ReferenceBases>> refCollection = p.apply(Create.of(refSource.getAllReferenceBases()));
        final PCollectionView<Map<String, ReferenceBases>> refView = refCollection.apply(View.asMap());

        p.apply(Create.<Void>of((Void) null))
                .setCoder(VoidCoder.of())
                .apply(ParDo.of(new DoFn<Void, Void>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        Map<String, ReferenceBases> ref = c.sideInput(refView);
                        Assert.assertNotNull(ref);
                        Assert.assertEquals(ref.size(), 4);
                        ReferenceBases bases = ref.get("2");
                        ReferenceBases subset = bases.getSubset(new SimpleInterval("2", 10001, 10010));
                        Assert.assertNotNull(subset);
                        Assert.assertEquals(subset.getBases().length, 10, "Wrong number of bases returned");
                        Assert.assertEquals(new String(subset.getBases()), "CGTATCCCAC", "Wrong bases returned");
                    }
                }).withSideInputs(refView));

        p.run();

    }

}
