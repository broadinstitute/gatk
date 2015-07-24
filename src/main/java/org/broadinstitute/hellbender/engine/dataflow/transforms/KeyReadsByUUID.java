package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.UUID;

/**
 * KeyReadsByUUID takes a PCollection of reads and produces one KV per read where the key is the UUID from the read.
 * This is necessary to use reads in GroupBy and CoGroupBy since GATKRead cannot be deterministically encoded.
 */
public class KeyReadsByUUID extends PTransform<PCollection<GATKRead>, PCollection<KV<UUID, GATKRead>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<UUID, GATKRead>> apply(PCollection<GATKRead> input) {
        return input.apply(ParDo.of(new DoFn<GATKRead, KV<UUID, GATKRead>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                c.output(KV.of(c.element().getUUID(), c.element()));
            }
        })).setName("KeyReadsByUUID");
    }
}
