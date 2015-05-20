package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.Read;

import java.util.List;

/**
 * Created by davidada on 5/15/15.
 */
public class KeyReadByShard extends PTransform<PCollection<Read>, PCollection<KV<SimpleInterval, Read>>> {
    @Override
    public PCollection<KV<SimpleInterval, Read>> apply(PCollection<Read> input) {
        return input.apply(ParDo.of(new DoFn<Read, KV<SimpleInterval, Read>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                List<SimpleInterval> intervals = DataflowUtils.getVariantShardsFromInterval(c.element());
                for (SimpleInterval inteval : intervals) {
                    c.output(KV.of(inteval, c.element()));
                }
            }
        }).named("KeyLocatableByShard"));
    }

}
