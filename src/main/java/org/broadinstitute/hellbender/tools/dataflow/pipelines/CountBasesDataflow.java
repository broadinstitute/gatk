package org.broadinstitute.hellbender.tools.dataflow.pipelines;


import com.google.common.collect.ImmutableList;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.dataflow.transforms.CountBasesDataflowTransform;

import java.util.logging.Logger;


@CommandLineProgramProperties(usage = "Count bases in dataflow with filters and transformers", usageShort = "count bases", programGroup = DataFlowProgramGroup.class)
public final class CountBasesDataflow extends DataflowReadsPipeline {
    @Override
    protected ImmutableList<ReadFilter> getReadFilters(){
        return ImmutableList.of(ReadFilterLibrary.WELLFORMED);
    }

    @Override
    protected PTransformSAM getTool() {
        return new CountBasesDataflowTransform();
    }

}



