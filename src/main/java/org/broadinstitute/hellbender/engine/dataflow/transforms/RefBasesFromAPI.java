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
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

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
 *
 * If a custom reference window function is being used to map each read to arbitrary reference bases
 * (and not just the bases that overlap each read), that function should be passed in via the
 * {@link ReferenceDataflowSource} parameter to {@link #getBasesForShard}.
 */
public class RefBasesFromAPI {
    public static PCollection<KV<ReferenceBases, Iterable<GATKRead>>> getBasesForShard(PCollection<KV<ReferenceShard, Iterable<GATKRead>>> reads,
                                                                                       ReferenceDataflowSource referenceDataflowSource) {
        PCollectionView<ReferenceDataflowSource> refView = reads.getPipeline().apply("apply create of referenceDataflowSource",
                Create.of(referenceDataflowSource)).apply("View ReferenceDataflowSource as singleton",View.<ReferenceDataflowSource>asSingleton());
        return reads.apply(ParDo.withSideInputs(refView).of(
                new DoFn<KV<ReferenceShard, Iterable<GATKRead>>, KV<ReferenceBases, Iterable<GATKRead>>>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(ProcessContext c) throws Exception {
                        final Iterable<GATKRead> reads = c.element().getValue();

                        // Apply the reference window function to each read to produce a set of intervals representing
                        // the desired reference bases for each read.
                        final List<SimpleInterval> readWindows = StreamSupport.stream(reads.spliterator(), false).map(read -> referenceDataflowSource.getReferenceWindowFunction().apply(read)).collect(Collectors.toList());

                        // Get a single interval spanning all the per-read reference windows.
                        SimpleInterval interval = SimpleInterval.getSpanningInterval(readWindows);
                        ReferenceBases bases = c.sideInput(refView).getReferenceBases(interval, c.getPipelineOptions());
                        c.output(KV.of(bases, reads));
                    }
                })).setName("RefBasesFromAPI_getBasesForShard");
    }
}
