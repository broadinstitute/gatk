package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;

public class ReadToInterval extends PTransform<PCollection<Read>, PCollection<SimpleInterval>> {

    @Override
    public PCollection<SimpleInterval> apply( PCollection<Read> input ) {
        return input.apply(ParDo.of(new DoFn<Read, SimpleInterval>() {
            @Override
            public void processElement( ProcessContext c ) throws Exception {
                c.output(new SimpleInterval(c.element()));
            }
        }).named("ReadToInterval"));
    }
}

