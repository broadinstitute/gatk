package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public final class RemoveDuplicateReadVariantPairsTest {

    @DataProvider(name = "dupedPairedReadsAndVariants")
    public Object[][] dupedPairedReadsAndVariants() {
        DataflowTestData testData = new DataflowTestData();
        List<KV<GATKRead, Variant>> dupes = testData.getKvReadVariant();

        List<KV<UUID, UUID>> kvUUIDUUID = Arrays.asList(
                KV.of(dupes.get(0).getKey().getUUID(), dupes.get(0).getValue().getUUID()),
                KV.of(dupes.get(1).getKey().getUUID(), dupes.get(1).getValue().getUUID()),
                KV.of(dupes.get(2).getKey().getUUID(), dupes.get(2).getValue().getUUID()),
                KV.of(dupes.get(3).getKey().getUUID(), dupes.get(3).getValue().getUUID()),
                KV.of(dupes.get(4).getKey().getUUID(), dupes.get(4).getValue().getUUID()));

        Iterable<UUID> uuids0 = Arrays.asList(dupes.get(0).getValue().getUUID(), dupes.get(1).getValue().getUUID());
        Iterable<UUID> uuids2 = Arrays.asList(dupes.get(2).getValue().getUUID());
        Iterable<UUID> uuids3 = Arrays.asList(dupes.get(3).getValue().getUUID());
        List<KV<UUID, Iterable<UUID>>> kvUUIDiUUID = Arrays.asList(
                KV.of(dupes.get(0).getKey().getUUID(), uuids0),
                KV.of(dupes.get(2).getKey().getUUID(), uuids2),
                KV.of(dupes.get(3).getKey().getUUID(), uuids3)
        );



        Iterable<Variant> variant01 = Arrays.asList(dupes.get(1).getValue(), dupes.get(0).getValue());
        Iterable<Variant> variant2 = Arrays.asList(dupes.get(2).getValue());
        Iterable<Variant> variant3 = Arrays.asList(dupes.get(3).getValue());
        List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant = Arrays.asList(
                KV.of(dupes.get(0).getKey().getUUID(), variant01),
                KV.of(dupes.get(2).getKey().getUUID(), variant2),
                KV.of(dupes.get(3).getKey().getUUID(), variant3));

        List<KV<GATKRead, Iterable<Variant>>> finalExpected = testData.getKvReadiVariant();

        return new Object[][]{
                {dupes, kvUUIDUUID, kvUUIDiUUID, kvUUIDiVariant, finalExpected},
        };
    }

    @Test(dataProvider = "dupedPairedReadsAndVariants")
    public void removeDupedReadVariantsTest(List<KV<GATKRead, Variant>> dupes, List<KV<UUID, UUID>> kvUUIDUUID,
                                            List<KV<UUID, Iterable<UUID>>> kvUUIDiUUID, List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant,
                                            List<KV<GATKRead, Iterable<Variant>>> finalExpected) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<GATKRead, Variant>> pKVs = DataflowTestUtils.PCollectionCreateAndVerify(p, dupes);

        PCollection<KV<UUID, UUID>> uuids = pKVs.apply(ParDo.of(new RemoveDuplicateReadVariantPairs.StripToUUIDs()));
        DataflowAssert.that(uuids).containsInAnyOrder(kvUUIDUUID);

        PCollection<KV<UUID, Iterable<UUID>>> readsIterableVariants = uuids.apply(RemoveDuplicates.<KV<UUID, UUID>>create()).apply(GroupByKey.<UUID, UUID>create());
        readsIterableVariants.apply(ParDo.of(new PrintUUIDs()));
        DataflowAssert.that(readsIterableVariants).containsInAnyOrder(kvUUIDiUUID);

        // Now add the variants back in.
        // Now join the stuff we care about back in.
        // Start by making KV<UUID, KV<UUID, Variant>> and joining that to KV<UUID, Iterable<UUID>>
        PCollection<KV<UUID, Iterable<Variant>>> matchedVariants = RemoveDuplicateReadVariantPairs.RemoveReadVariantDupesUtility.addBackVariants(pKVs, readsIterableVariants);
        DataflowAssert.that(matchedVariants).containsInAnyOrder(kvUUIDiVariant);

        PCollection<KV<GATKRead, Iterable<Variant>>> finalResult = RemoveDuplicateReadVariantPairs.RemoveReadVariantDupesUtility.addBackReads(pKVs, matchedVariants);
        DataflowAssert.that(finalResult).containsInAnyOrder(finalExpected);

        p.run();
    }

    @Test(dataProvider = "dupedPairedReadsAndVariants")
    public void fullRemoveDupesTest(List<KV<GATKRead, Variant>> dupes, List<KV<UUID, UUID>> kvUUIDUUID,
                                            List<KV<UUID, Iterable<UUID>>> kvUUIDiUUID, List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant,
                                            List<KV<GATKRead, Iterable<Variant>>> finalExpected) {
        Pipeline p = TestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<GATKRead, Variant>> pKVs = DataflowTestUtils.PCollectionCreateAndVerify(p, dupes);

        PCollection<KV<GATKRead, Iterable<Variant>>> result = pKVs.apply(new RemoveDuplicateReadVariantPairs());
        DataflowAssert.that(result).containsInAnyOrder(finalExpected);

        p.run();
    }

    static class PrintUUIDs extends DoFn<KV<UUID, Iterable<UUID>>, Integer> {
        private static final long serialVersionUID = 1L;
        @Override
        public void processElement(ProcessContext c) throws Exception {
            KV<UUID, Iterable<UUID>> i = c.element();
            System.out.print("--- " + i.getKey() + ": ");
            for (UUID uuid : i.getValue()) {
                System.out.print(uuid + ",");
            }
            System.out.println(" ---");
        }
    }

}