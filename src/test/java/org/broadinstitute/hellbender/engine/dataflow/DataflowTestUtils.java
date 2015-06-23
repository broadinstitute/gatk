package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.DataflowAssert;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.Lists;
import org.testng.Assert;

public class DataflowTestUtils {
    public static <T> PCollection<T> PCollectionCreateAndVerify(Pipeline p, Iterable<T> list) {
        Iterable<T> copy = Lists.newArrayList(list.iterator());
        Assert.assertEquals(list, copy);
        PCollection<T> pCollection = p.apply(Create.of(list));
        DataflowAssert.that(pCollection).containsInAnyOrder(copy);
        return pCollection;
    }
}
