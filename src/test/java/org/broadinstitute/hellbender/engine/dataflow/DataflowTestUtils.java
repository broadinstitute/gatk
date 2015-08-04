package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.Coder;
import com.google.cloud.dataflow.sdk.testing.CoderProperties;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.common.collect.Lists;
import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;

import java.util.List;

public class DataflowTestUtils {
    /**
     * PCollectionCreateAndVerify asserts that the contents of the collection can be created, then returns the resulting
     * PCollection.
     *
     * @param p     to use with Create.of()
     * @param list  the items to put into the PCollection
     * @param coder the coder for the items, new GATKReadCoder.
     * @param <T>   the item type, e.g., GATKRead
     * @return a PCollection containing the items from list.
     */
    public static <T> PCollection<T> pCollectionCreateAndVerify(Pipeline p, Iterable<T> list, Coder<T> coder) {
        for (T value : list) {
            try {
                CoderProperties.coderDecodeEncodeEqual(coder, value);
            } catch (Exception e) {
                throw new GATKException("DecodeEncodeEqual are not equal for: " + value.getClass().getSimpleName());
            }
        }
        PCollection<T> pCollection = p.apply(Create.of(list).withCoder(coder));
        DataflowAssert.that(pCollection).containsInAnyOrder(list); // Remove me when DavidR's fix is in.
        return pCollection;
    }

    static class BInA<T, U> extends DoFn<KV<T, Iterable<U>>, Void> {
        private static final long serialVersionUID = 1L;
        PCollectionView<List<KV<T, Iterable<U>>>> view;
        public BInA(PCollectionView<List<KV<T, Iterable<U>>>> view) {
            super();
            this.view = view;
        }

        @Override
        public void processElement(ProcessContext c) throws Exception {
            List<KV<T, Iterable<U>>> kvs2 = c.sideInput(view);
            KV<T, Iterable<U>> kv1 = c.element();
            T k1 = kv1.getKey();
            Iterable<U> v1 = kv1.getValue();
            boolean foundMatch = false;
            for (KV<T, Iterable<U>> kv2 : kvs2) {
                if (k1.equals(kv2.getKey())) {
                    // Found the key, now check the values (in any permutation).
                    if (CollectionUtils.isEqualCollection(Lists.newArrayList(v1), Lists.newArrayList(kv2.getValue()))) {
                        // Found a match!
                        foundMatch = true;
                    }
                }
            }
            Assert.assertTrue(foundMatch, "Unable to find a match for " + kv1.toString() + " in p2");
        }
    }

    public static <T, U> void keyIterableValueMatcher(PCollection<KV<T, Iterable<U>>> p1, PCollection<KV<T, Iterable<U>>> p2) {
        PCollectionView<List<KV<T, Iterable<U>>>> p2View = p2.apply(View.asList());
        p1.apply(ParDo.withSideInputs(p2View).of(new BInA<>(p2View)));

        PCollectionView<List<KV<T, Iterable<U>>>> p1View = p2.apply("p1View.View",View.asList());
        p2.apply("p2.ParDo.BInA",ParDo.withSideInputs(p1View).of(new BInA<>(p1View)));
    }

    static class BInAwReadContextData<T> extends DoFn<KV<T, ReadContextData>, Void> {
        private static final long serialVersionUID = 1L;
        PCollectionView<List<KV<T, ReadContextData>>> view;
        public BInAwReadContextData(PCollectionView<List<KV<T, ReadContextData>>> view) {
            super();
            this.view = view;
        }

        @Override
        public void processElement(ProcessContext c) throws Exception {
            List<KV<T, ReadContextData>> kvs2 = c.sideInput(view);
            KV<T, ReadContextData> kv1 = c.element();
            T k1 = kv1.getKey();
            ReadContextData v1 = kv1.getValue();
            boolean foundMatch = false;
            for (KV<T, ReadContextData> kv2 : kvs2) {
                if (k1.equals(kv2.getKey())) {
                    // Found the key, now check the ReadContextData (Ref bases with Variants in any permutation).
                    if (CollectionUtils.isEqualCollection(Lists.newArrayList(v1.getOverlappingVariants()),
                            Lists.newArrayList(kv2.getValue().getOverlappingVariants())) &&
                            v1.getOverlappingReferenceBases().equals(kv2.getValue().getOverlappingReferenceBases())) {
                        // Found a match!
                        foundMatch = true;
                    }
                }
            }
            Assert.assertTrue(foundMatch, "Unable to find a match for " + kv1.toString() + " in p2");
        }
    }

    public static <T> void keyReadContextDataMatcher(PCollection<KV<T, ReadContextData>> p1, PCollection<KV<T, ReadContextData>> p2) {
        PCollectionView<List<KV<T, ReadContextData>>> p2View = p2.apply(View.asList());
        p1.apply(ParDo.withSideInputs(p2View).of(new BInAwReadContextData<>(p2View)));

        PCollectionView<List<KV<T, ReadContextData>>> p1View = p2.apply(View.asList());
        p2.apply(ParDo.withSideInputs(p1View).of(new BInAwReadContextData<>(p1View)));
    }
}
