package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * RefBasesFromAPI queries the Google Genomics API for reference bases that span all of the reads on each shard.
 *
 * |---------------- shard 0 ----------------|---------------- shard 1 ----------------|
 * |-- read a --|   |---- read b ---|           |--------------- read c --------------|
 *
 * results in,
 * |---------- ref bases 1 ---------|           |---------- ref bases 2 --------------|
 * |-- read a --|   |---- read b ---|           |--------------- read c --------------|
 *
 * KV<ref bases 1, [read a, read b]>
 * KV<ref bases 2, [read c]>
 */
public class RefBasesFromAPI {
    public static PCollection<KV<ReferenceBases, Iterable<GATKRead>>> getBasesForShard(PCollection<KV<ReferenceShard, Iterable<GATKRead>>> reads,
                                                                                       RefAPIMetadata refAPIMetadata) {
        PCollectionView<RefAPIMetadata> dataView = reads.getPipeline().apply(Create.of(refAPIMetadata)).apply(View.<RefAPIMetadata>asSingleton());
        return reads.apply(ParDo.withSideInputs(dataView).of(
                new DoFn<KV<ReferenceShard, Iterable<GATKRead>>, KV<ReferenceBases, Iterable<GATKRead>>>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        final Iterable<GATKRead> reads = c.element().getValue();
                        SimpleInterval interval = SimpleInterval.getSpanningInterval(reads);
                        RefAPISource refAPISource = RefAPISource.getRefAPISource();

                        ReferenceBases bases = refAPISource.getReferenceBases(c.getPipelineOptions(), c.sideInput(dataView), interval);
                        c.output(KV.of(bases, reads));
                    }
                })).setName("RefBasesFromAPI_getBasesForShard");
    }
}
