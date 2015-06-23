package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;

/**
 * transformExample is a simple transform that takes Integers and produces a KV pair of ReferenceShard and the same
 * Integer. The purpose is to show the minimum pieces required to make a PTransform.
 */
public class transformExample extends PTransform<PCollection<Integer>, PCollection<KV<ReferenceShard, Integer>>> {
    private static final long serialVersionUID = 1L;
    @Override
    public PCollection<KV<ReferenceShard, Integer>> apply(PCollection<Integer> input) {
            return input.apply(ParDo.of(new DoFn<Integer, KV<ReferenceShard, Integer>>() {
                private static final long serialVersionUID = 1L;
                @Override
                public void processElement(ProcessContext c) throws Exception {
                    Integer i = c.element();
                    c.output(KV.of(new ReferenceShard(i.intValue(), i.toString()), i));
                }
            }).named("transformExample"));
    }
}

