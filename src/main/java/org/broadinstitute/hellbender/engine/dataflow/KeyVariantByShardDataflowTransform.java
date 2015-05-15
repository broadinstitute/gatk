package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;

/**
 * Created by davidada on 5/15/15.
 */
public class KeyVariantByShardDataflowTransform extends PTransform<PCollection<Variant>, PCollection<KV<SimpleInterval, Variant>>> {
    @Override
    public PCollection<KV<SimpleInterval, Variant>> apply(PCollection<Variant> input) {
        return input.apply(ParDo.of(new DoFn<Variant, KV<SimpleInterval, Variant>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                List<SimpleInterval> intervals = DataflowUtils.getShardsFromInterval(c.element());
                for (SimpleInterval inteval : intervals) {
                    c.output(KV.of(inteval, c.element()));
                }
            }
        }).named("KeyVariantByShard"));
    }

}
