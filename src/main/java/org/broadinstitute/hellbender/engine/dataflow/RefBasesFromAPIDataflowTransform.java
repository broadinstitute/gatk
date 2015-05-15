package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.Iterator;

/**
 * Created by davidada on 5/15/15.
 */
public class RefBasesFromAPIDataflowTransform extends PTransform<PCollection<KV<ReferenceShard, Read>>, PCollection<KV<ReferenceBases, Iterable<Read>>>> {
    @Override
    public PCollection<KV<ReferenceBases, Iterable<Read>>> apply(PCollection<KV<ReferenceShard, Read>> input) {

        PCollection<KV<ReferenceShard, Iterable<Read>>> keyed = input.apply(GroupByKey.<ReferenceShard, Read>create());

        return keyed.apply(ParDo.of(new DoFn<KV<ReferenceShard, Iterable<Read>>, KV<ReferenceBases, Iterable<Read>>>() {
            @Override
            public void processElement( ProcessContext c ) throws Exception {
                final ReferenceShard shard = c.element().getKey();
                final Iterable<Read> reads = c.element().getValue();
                int min = 1000000000; // ????
                int max = 0;
                for (Read r : reads) {
                    if (r.getStart() < min) {
                        min = r.getStart();
                    }
                    if (r.getEnd() < max) {
                        max = r.getEnd();
                    }
                }
                ReferenceSource source = new ReferenceSource("ref_name", c.getPipelineOptions());
                ReferenceBases bases = source.getReferenceBases(new SimpleInterval(shard.getContig(), min, max));
                c.output(KV.of(bases, reads));
            }
        }));
    }
}
