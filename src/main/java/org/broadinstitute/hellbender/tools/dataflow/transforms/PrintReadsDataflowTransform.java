package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowSAMFn;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;


/**
 * Print a string representation of reads in the PCollection<String>
 */
public final class PrintReadsDataflowTransform extends PTransformSAM<String> {
        private static final long serialVersionUID = 1l;

        @Override
        public PCollection<String> apply(final PCollection<Read> reads) {
                return reads.apply(ParDo.of(new DataFlowSAMFn<String>(getHeader()) {
                        private static final long serialVersionUID = 1l;

                        @Override
                        protected void apply(final SAMRecord read) {
                                output(read.getSAMString());
                        }
                }));
        }
}
