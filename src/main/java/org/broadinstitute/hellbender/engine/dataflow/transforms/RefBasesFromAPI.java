package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * Created by davidada on 5/15/15.
 */
public class RefBasesFromAPI extends PTransform<PCollection<KV<ReferenceShard, Read>>, PCollection<KV<ReferenceBases, Iterable<Read>>>> {
    private final String refName;

    public RefBasesFromAPI( String refName ) {
        this.refName = refName;
    }

    @Override
    public PCollection<KV<ReferenceBases, Iterable<Read>>> apply(PCollection<KV<ReferenceShard, Read>> input) {
        PCollection<KV<ReferenceShard, Iterable<Read>>> keyed = input.apply(GroupByKey.<ReferenceShard, Read>create());
        PCollectionView<String> refView = input.getPipeline().apply(Create.of(refName)).apply(View.<String>asSingleton());

        return keyed.apply(ParDo.withSideInputs(refView).of(new DoFn<KV<ReferenceShard, Iterable<Read>>, KV<ReferenceBases, Iterable<Read>>>() {
            @Override
            public void processElement(ProcessContext c) throws Exception {
                final ReferenceShard shard = c.element().getKey();
                final Iterable<Read> reads = c.element().getValue();
                int min = Integer.MAX_VALUE;
                int max = 0;
                for (Read r : reads) {
                    if (r.getStart() < min) {
                        min = r.getStart();
                    }
                    if (r.getEnd() < max) {
                        max = r.getEnd();
                    }
                }
                ReferenceSource source = new ReferenceSource(c.sideInput(refView), c.getPipelineOptions());
                ReferenceBases bases = source.getReferenceBases(new SimpleInterval(shard.getContig(), min, max));
                c.output(KV.of(bases, reads));
            }
        }));
    }
}
