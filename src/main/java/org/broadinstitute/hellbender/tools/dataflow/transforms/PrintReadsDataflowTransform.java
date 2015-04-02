package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowReadFn;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.utils.read.GATKRead;


/**
 * Print a string representation of reads in the PCollection<String>
 */
public final class PrintReadsDataflowTransform extends PTransformSAM<String> {
        private static final long serialVersionUID = 1l;

        @Override
        public PCollection<String> apply(final PCollection<GATKRead> reads) {
                return reads.apply(ParDo.of(new DataFlowReadFn<String>(getHeader()) {
                        private static final long serialVersionUID = 1l;

                        @Override
                        protected void apply(final GATKRead read) {
                                // TODO: write a utility that can produce a SAM string for a GATK read without conversion to SAMRecord
                                output(read.convertToSAMRecord(getHeader()).getSAMString());
                        }
                }));
        }
}
