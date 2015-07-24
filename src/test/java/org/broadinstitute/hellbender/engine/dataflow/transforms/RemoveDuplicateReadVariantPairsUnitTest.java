package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.IterableCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.GATKTestPipeline;
import org.broadinstitute.hellbender.engine.dataflow.ReadsPreprocessingPipelineTestData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTestUtils;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.UUIDCoder;
import org.broadinstitute.hellbender.engine.dataflow.coders.VariantCoder;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import java.util.*;

public final class RemoveDuplicateReadVariantPairsUnitTest extends BaseTest {

    @DataProvider(name = "dupedPairedReadsAndVariants")
    public Object[][] dupedPairedReadsAndVariants() {
        Object[][] data = new Object[2][];
        List<Class<?>> classes = Arrays.asList(Read.class, SAMRecord.class);
        for (int i = 0; i < classes.size(); ++i) {
            Class<?> c = classes.get(i);
            ReadsPreprocessingPipelineTestData testData = new ReadsPreprocessingPipelineTestData(c);
            List<KV<GATKRead, Variant>> dupes = testData.getKvReadVariant();

            List<KV<UUID, UUID>> kvUUIDUUID = Arrays.asList(
                    KV.of(dupes.get(0).getKey().getUUID(), dupes.get(0).getValue().getUUID()),
                    KV.of(dupes.get(1).getKey().getUUID(), dupes.get(1).getValue().getUUID()),
                    KV.of(dupes.get(2).getKey().getUUID(), dupes.get(2).getValue().getUUID()),
                    KV.of(dupes.get(3).getKey().getUUID(), dupes.get(3).getValue().getUUID()),
                    KV.of(dupes.get(4).getKey().getUUID(), dupes.get(4).getValue().getUUID()),
                    KV.of(dupes.get(5).getKey().getUUID(), dupes.get(5).getValue().getUUID())
            );

            Iterable<UUID> uuids0 = Arrays.asList(dupes.get(1).getValue().getUUID(), dupes.get(0).getValue().getUUID());
            Iterable<UUID> uuids2 = Arrays.asList(dupes.get(2).getValue().getUUID());
            Iterable<UUID> uuids3 = Arrays.asList(dupes.get(3).getValue().getUUID());
            Iterable<UUID> uuids5 = Arrays.asList(dupes.get(5).getValue().getUUID());
            List<KV<UUID, Iterable<UUID>>> kvUUIDiUUID = Arrays.asList(
                    KV.of(dupes.get(0).getKey().getUUID(), uuids0),
                    KV.of(dupes.get(2).getKey().getUUID(), uuids2),
                    KV.of(dupes.get(3).getKey().getUUID(), uuids3),
                    KV.of(dupes.get(5).getKey().getUUID(), uuids5)
            );

            Iterable<Variant> variant01 = Arrays.asList(dupes.get(1).getValue(), dupes.get(0).getValue());
            Iterable<Variant> variant2 = Arrays.asList(dupes.get(2).getValue());
            Iterable<Variant> variant3 = Arrays.asList(dupes.get(3).getValue());
            Iterable<Variant> variant5 = Arrays.asList(dupes.get(5).getValue());

            List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant = Arrays.asList(
                    KV.of(dupes.get(0).getKey().getUUID(), variant01),
                    KV.of(dupes.get(2).getKey().getUUID(), variant2),
                    KV.of(dupes.get(3).getKey().getUUID(), variant3),
                    KV.of(dupes.get(5).getKey().getUUID(), variant5)
            );

            List<KV<GATKRead, Iterable<Variant>>> finalExpected = testData.getKvReadiVariant();
            data[i] = new Object[]{dupes, kvUUIDUUID, kvUUIDiUUID, kvUUIDiVariant, finalExpected};
        }
        return data;
    }

    @Test(dataProvider = "dupedPairedReadsAndVariants")
    public void removeDupedReadVariantsTest(List<KV<GATKRead, Variant>> dupes, List<KV<UUID, UUID>> kvUUIDUUID,
                                            List<KV<UUID, Iterable<UUID>>> kvUUIDiUUID, List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant,
                                            List<KV<GATKRead, Iterable<Variant>>> finalExpected) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<GATKRead, Variant>> pKVs = DataflowTestUtils.pCollectionCreateAndVerify(p, dupes,
                KvCoder.of(new GATKReadCoder(), new VariantCoder()));

        PCollection<KV<UUID, UUID>> uuids = pKVs.apply(ParDo.of(new RemoveDuplicateReadVariantPairs.StripToUUIDs()));
        DataflowAssert.that(uuids).containsInAnyOrder(kvUUIDUUID);

        PCollection<KV<UUID, Iterable<UUID>>> readsIterableVariants = uuids.apply(RemoveDuplicates.<KV<UUID, UUID>>create()).apply(GroupByKey.<UUID, UUID>create());
        PCollection<KV<UUID, Iterable<UUID>>> tempExpected = p.apply(Create.of(kvUUIDiUUID)).setCoder(KvCoder.of(UUIDCoder.CODER, IterableCoder.of(UUIDCoder.CODER)));
        DataflowTestUtils.keyIterableValueMatcher(readsIterableVariants, tempExpected);

        // Now add the variants back in.
        // Now join the stuff we care about back in.
        // Start by making KV<UUID, KV<UUID, Variant>> and joining that to KV<UUID, Iterable<UUID>>
        PCollection<KV<UUID, Iterable<Variant>>> matchedVariants = RemoveDuplicateReadVariantPairs.RemoveReadVariantDupesUtility.addBackVariants(pKVs, readsIterableVariants);
        PCollection<KV<UUID, Iterable<Variant>>> pkvUUIDiVariant = p.apply(Create.of(kvUUIDiVariant)).setCoder(KvCoder.of(UUIDCoder.CODER, IterableCoder.of(new VariantCoder())));
        DataflowTestUtils.keyIterableValueMatcher(matchedVariants, pkvUUIDiVariant);

        PCollection<KV<GATKRead, Iterable<Variant>>> finalResult = RemoveDuplicateReadVariantPairs.RemoveReadVariantDupesUtility.addBackReads(pKVs, matchedVariants);
        PCollection<KV<GATKRead, Iterable<Variant>>> pFinalExpected = p.apply(Create.of(finalExpected)).setCoder(KvCoder.of(new GATKReadCoder(), IterableCoder.of(new VariantCoder())));
        DataflowTestUtils.keyIterableValueMatcher(finalResult, pFinalExpected);
        p.run();
    }

    @Test(dataProvider = "dupedPairedReadsAndVariants")
    public void fullRemoveDupesTest(List<KV<GATKRead, Variant>> dupes, List<KV<UUID, UUID>> kvUUIDUUID,
                                            List<KV<UUID, Iterable<UUID>>> kvUUIDiUUID, List<KV<UUID, Iterable<Variant>>> kvUUIDiVariant,
                                            List<KV<GATKRead, Iterable<Variant>>> finalExpected) {
        Pipeline p = GATKTestPipeline.create();
        DataflowUtils.registerGATKCoders(p);

        PCollection<KV<GATKRead, Variant>> pKVs = DataflowTestUtils.pCollectionCreateAndVerify(p, dupes,
                KvCoder.of(new GATKReadCoder(), new VariantCoder()));

        PCollection<KV<GATKRead, Iterable<Variant>>> result = pKVs.apply(new RemoveDuplicateReadVariantPairs());
        PCollection<KV<GATKRead, Iterable<Variant>>> pFinalExpected = p.apply(Create.of(finalExpected)).setCoder(KvCoder.of(new GATKReadCoder(), IterableCoder.of(new VariantCoder())));
        DataflowTestUtils.keyIterableValueMatcher(result, pFinalExpected);

        p.run();
    }
}