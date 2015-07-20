package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.tools.dataflow.transforms.PrintReadsDataflowTransform;

@CommandLineProgramProperties(
        summary = "Prints read from a bam file to a new text file.",
        oneLineSummary = "Print reads",
        programGroup = DataFlowProgramGroup.class
)
public final class PrintReadsDataflow extends DataflowReadsPipeline{
        private static final long serialVersionUID = 1l;

        @Override
        protected PTransformSAM<String> getTool() {
                return new PrintReadsDataflowTransform();
        }
}
