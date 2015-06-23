package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.List;

public final class ReaderCoderTest {

    @DataProvider(name = "reads")
    public Object[][] uuids(){

        List<GATKRead> googleReads;// = Lists.newArrayList(makeGoogleRead(1, 10, 1), makeGoogleRead(20, 40, 2));
        List<GATKRead> samReads = Lists.newArrayList(makeSamRead(1, 10, 1), makeSamRead(20, 40, 2));
        /*
        List<UUID> oneHundredUuids = Lists.newArrayList();
        for (int i = 0; i < 100; ++i) {
            oneHundredUuids.add(UUID.randomUUID());
        }*/

        return new Object[][]{
                {samReads/*, googleReads*/},
        };
    }

    /*
    private GATKRead makeGoogleRead(int start, int length, int i) {
        return ArtificialReadUtils.createRandomGoogleRead(start, length, i);
    }*/

    public GATKRead makeSamRead(int start, int length, int i) {
        return ArtificialReadUtils.createRandomRead(start, length, i);
    }

    @Test(dataProvider = "reads")
    public void createUuidsTest(List<GATKRead> samReads/*, List<GATKRead> googleReads*/) {
        // The simplest way to figure out if a class is coded correctly is to create a PCollection
        // of that type and see if matches the List version.
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);
        //p.getCoderRegistry().registerCoder(GATKRead.class, ReadCoder.CODER);
        //p.getCoderRegistry().registerCoder(GoogleGenomicsReadToGATKReadAdapter.class, GoogleGenomicsReadToGATKReadAdapter.CODER);
        //p.getCoderRegistry().registerCoder(SAMRecordToGATKReadAdapter.class, SerializableCoder.of(SAMRecordToGATKReadAdapter.class));

        PCollection<GATKRead> pShards = p.apply(Create.of(samReads));

        List<GATKRead> sameReads = Lists.newArrayList();
        Assert.assertTrue(sameReads.addAll(samReads));
        DataflowAssert.that(pShards).containsInAnyOrder(samReads);
        p.run();
    }

    /*
    @Test(dataProvider = "reads")
    public void uniqueUuidsTest(List<MutableRead> twoReads) {
        // Note that this is a probabilistic test. It's possible that this test could fail by chance,
        // but the test is designed to fail with fewer than a one in a million attempts (easily).
        Pipeline p = TestPipeline.create();
        p.getCoderRegistry().registerCoder(UUID.class, UuidCoder.CODER);

        PCollection<UUID> pShards = p.apply(Create.of(uuids));
        // createUuidsTest makes sure the UUIDs are coded correctly, we assume they are here.

        // Count the number of unique UUIDs and make sure that count matches the total number of UUIDs.
        PCollection<Long> counts = pShards.apply(RemoveDuplicates.<UUID>create()).apply(Count.<UUID>globally());
        List<Long> expectedCounts = Lists.newArrayList((long) uuids.size());
        DataflowAssert.that(counts).containsInAnyOrder(expectedCounts);
        p.run();
    }*/
}
