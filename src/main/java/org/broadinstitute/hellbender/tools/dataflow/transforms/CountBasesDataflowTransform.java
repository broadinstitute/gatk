package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Sum;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowReadFn;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Count the number of bases in a PCollection<GATKRead>
 */
public final class CountBasesDataflowTransform extends PTransformSAM<Long> {
    private static final long serialVersionUID = 1l;

    @Override
    public PCollection<Long> apply(final PCollection<GATKRead> reads) {

        return reads.apply(ParDo.of(new DataFlowReadFn<Long>(getHeader()) {
            private static final long serialVersionUID = 1l;

            @Override
            protected void apply(final GATKRead read) {
                final long bases = read.getLength();
                output(bases);
            }
        }))
        .apply(Sum.longsGlobally());
    }

}
