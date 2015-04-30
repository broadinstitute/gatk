package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;

public class KeyReadsByIntervalDataflowTransform extends PTransform<PCollection<Read>, PCollection<KV<SimpleInterval, Read>>> {

    @Override
    public PCollection<KV<SimpleInterval, Read>> apply( PCollection<Read> input ) {
        return input.apply(ParDo.of(new DoFn<Read, KV<SimpleInterval, Read>>() {
            @Override
            public void processElement( ProcessContext c ) throws Exception {
                final Read read = c.element();
                c.output(KV.of(new SimpleInterval(read), read));
            }
        }).named("KeyReadsByInterval"));
    }
}
