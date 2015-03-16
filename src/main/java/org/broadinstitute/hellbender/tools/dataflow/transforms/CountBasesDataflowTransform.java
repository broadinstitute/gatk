package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Sum;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowSAMFn;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;

/**
 * Count the number of bases in a PCollection<Read>
 */
public final class CountBasesDataflowTransform extends PTransformSAM<Long> {
    @Override
    public PCollection<Long> apply(PCollection<Read> reads) {

        return reads.apply(ParDo.of(new DataFlowSAMFn<Long>(getHeaderString()) {
            @Override
            protected void apply(SAMRecord read) {
                Long bases = (long) read.getReadBases().length;
                output(bases);
            }
        }))
                .apply(Sum.longsGlobally());
    }

}
