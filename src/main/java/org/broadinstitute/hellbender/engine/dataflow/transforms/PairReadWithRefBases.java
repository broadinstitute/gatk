package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

/**
 * PairReadWithRefBases pairs the ReferenceBases covering the reads with each read and produces KV with read as
 * the key. The reference bases were obtained from the Google Genomics API and are guaranteed to span all reads that
 * are paired. See RefBasesFromAPI.getBasesForShard for details on how the ref bases are collected.
 *
 * |---------- ref bases 1 ---------|      |-------------------- ref bases 2 --------------------|
 * |-- read a --|   |---- read b ---|         |--------------- read c --------------|
 * results in,
 *
 * |-- read a --|   |---- read b ---|         |--------------- read c --------------|
 * |ref bases 1a|   |-ref bases 1b -|         |----------- ref bases 2c ------------|
 *
 * KV<read a, ref bases 1a>
 * KV<read b, ref bases 1b>
 * KV<read c, ref bases 2c>
 * Precondition: each refbases must contain ALL of the bases for all of the GATKReads it's paired with.
 *
 * Also, it's possible (and fine) for a read to span two ref bases. This happens when a read ends after the start of
 * first read on the next shard
 *
 * |------- shard 1 -------|------- shard 2 -------|
 * [---------- read a ---------{====]---- read 2 --------}
 *
 * |---------- ref bases 1 ---------|
 *                             |------ ref bases 2 ------|
 *
 * If a custom reference window function is being used to map each read to arbitrary reference bases
 * (and not just the bases that overlap each read), that function should be passed in at construction.
 */

public class PairReadWithRefBases extends PTransform<PCollection<KV<ReferenceBases, Iterable<GATKRead>>>, PCollection<KV<GATKRead, ReferenceBases>>> {
    private static final long serialVersionUID = 1L;

    private final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction;

    public PairReadWithRefBases() {
        this(RefWindowFunctions.IDENTITY_FUNCTION);
    }

    /**
     * @param referenceWindowFunction custom reference window function used to map each read to arbitrary reference bases
     */
    public PairReadWithRefBases( final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction ) {
        this.referenceWindowFunction = referenceWindowFunction;
    }

    @Override
    public PCollection<KV<GATKRead, ReferenceBases>> apply(PCollection<KV<ReferenceBases, Iterable<GATKRead>>> input) {
        return input.apply(ParDo.of(new DoFn<KV<ReferenceBases, Iterable<GATKRead>>, KV<GATKRead, ReferenceBases>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                // Each element of the PCollection is a set of reads keyed by a reference shard
                // The shard MUST have all of the reference bases for ALL of the reads. If not
                // it's an error.
                final ReferenceBases shard = c.element().getKey();
                final Iterable<GATKRead> reads = c.element().getValue();
                // For every read, find the subset of the reference that matches it according to our referenceWindowFunction
                for (GATKRead r : reads) {
                    final ReferenceBases subset = shard.getSubset(referenceWindowFunction.apply(r));
                    c.output(KV.of(r, subset));
                }
            }
        })).setName("GroupReadWithRefBases");
    }
}
