package org.broadinstitute.hellbender.tools.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Sum;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowSAMFn;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.tools.dataflow.pipelines.DataflowReadsPipeline;

@CommandLineProgramProperties(
        usage = "Prints reads from the input to the output.",
        usageShort = "Print reads",
        programGroup = DataFlowProgramGroup.class
)
public class PrintReadsDataflow extends DataflowReadsPipeline{

        @Override
        protected PTransformSAM getTool() {
                return new PrintReads();
        }


        public static class PrintReads extends PTransformSAM{
                @Override
                public PCollection<String> apply(PCollection<Read> reads) {
                        return reads.apply(ParDo.of(new DataFlowSAMFn<String>(getHeaderString()) {
                                @Override
                                protected void apply(SAMRecord read) {
                                        output(read.getSAMString());
                                }
                        }));
                }
        }
}
