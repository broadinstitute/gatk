package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.GroupByKey;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.Read;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * Created by davidada on 5/15/15.
 */
public class GroupReadWithRefBases extends PTransform<PCollection<KV<ReferenceBases, Iterable<Read>>>, PCollection<KV<Read, ReferenceBases>>> {
        @Override
        public PCollection<KV<Read, ReferenceBases>> apply(PCollection<KV<ReferenceBases, Iterable<Read>>> input) {
            return input.apply(ParDo.of(new DoFn<KV<ReferenceBases, Iterable<Read>>, KV<Read, ReferenceBases>>() {
                @Override
                public void processElement(ProcessContext c) throws Exception {
                    // Each element of the PCollection is a set of reads keyed by a reference shard
                    // The shard MUST have all of the reference bases for ALL of the reads. If not
                    // it's an error.
                    final ReferenceBases shard = c.element().getKey();
                    final Iterable<Read> reads = c.element().getValue();
                    // For every read, find the subset of the reference that matches it.
                    for (Read r : reads) {
                        final ReferenceBases subset = shard.getSubset(new SimpleInterval(r));
                        c.output(KV.of(r, subset));
                    }
                }
            }));
        }
    }
