package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.tools.dataflow.pipelines.DataflowReadsPipeline;
import org.broadinstitute.hellbender.tools.dataflow.transforms.CountReadsDataflowTransform;

@CommandLineProgramProperties(usage= "count reads using dataflow", usageShort = "count reads", programGroup = DataFlowProgramGroup.class)
public final class CountReadsDataflow extends DataflowReadsPipeline{

    @Override
    protected PTransformSAM<Long> getTool() {
        return new CountReadsDataflowTransform();
    }

}
