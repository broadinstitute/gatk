package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.Serializable;
import java.util.*;


public class GATKReadCoderUnitTest extends BaseTest implements Serializable {
    private static final long serialVersionUID = 1l;

    private GATKRead makeGoogleRead( final int uuid, final String name, final int start, final int length ) {
        return ArtificialReadUtils.createGoogleBackedReadWithUUID(new UUID(0, uuid), name, start, length);
    }

    private List<GATKRead> makeGoogleReads() {
        return Arrays.asList(makeGoogleRead(1, "google1", 1, 10), makeGoogleRead(2, "google2", 20, 40));
    }

    private GATKRead makeSamRead( final int uuid, final String name, final int start, final int length ) {
        return ArtificialReadUtils.createSamBackedReadWithUUID(new UUID(0, uuid), name, start, length);
    }

    private List<GATKRead> makeSamReads() {
        return Arrays.asList(makeSamRead(1, "sam1", 1, 10), makeSamRead(2, "sam2", 20, 40));
    }

    @DataProvider(name = "reads")
    public Object[][] makeReads() {
        final List<GATKRead> mixedReads = new ArrayList<>();
        mixedReads.addAll(makeGoogleReads());
        mixedReads.addAll(makeSamReads());

        return new Object[][] {
                { makeGoogleReads() },
                { makeSamReads() },
                { mixedReads }
        };
    }

    @Test(dataProvider = "reads")
    public void testGATKReadCoding( final List<GATKRead> reads ) {
        // The simplest way to figure out if a class is coded correctly is to create a PCollection
        // of that type and see if it matches the List version.
        final Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        // Need to explicitly set the coder to GATKReadCoder, otherwise Create fails to infer
        // a coder properly in the case where the List contains a mix of different GATKRead implementations.
        final PCollection<GATKRead> dataflowReads = p.apply(Create.of(reads)).setCoder(new GATKReadCoder());
        DataflowAssert.that(dataflowReads).containsInAnyOrder(reads);

        final PCollection<GATKRead> dataflowReadsAfterTransform = dataflowReads.apply(ParDo.of(new DoFn<GATKRead, GATKRead>() {
            private static final long serialVersionUID = 1l;

            @Override
            public void processElement( ProcessContext c ) throws Exception {
                c.output(c.element());
            }
        })).setCoder(new GATKReadCoder());
        DataflowAssert.that(dataflowReadsAfterTransform).containsInAnyOrder(reads);

        p.run();
    }
}
