package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * Created by davidada on 5/15/15.
 */
public class RefBasesToReadsDataflowTransform extends PTransform<PCollection<KV<ReferenceBases, Iterable<Read>>>, PCollection<KV<Read, ReferenceBases>>> {
        @Override
        public PCollection<KV<Read, ReferenceBases>> apply(PCollection<KV<ReferenceBases, Iterable<Read>>> input) {
            return input.apply(ParDo.of(new DoFn<KV<ReferenceBases, Iterable<Read>>, KV<Read, ReferenceBases>>() {
                @Override
                public void processElement(ProcessContext c) throws Exception {
                    final ReferenceBases shard = c.element().getKey();
                    final Iterable<Read> reads = c.element().getValue();
                    for (Read r : reads) {

                    }
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
