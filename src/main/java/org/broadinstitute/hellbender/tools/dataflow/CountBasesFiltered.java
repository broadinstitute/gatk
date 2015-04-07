package org.broadinstitute.hellbender.tools.dataflow;


import com.google.common.collect.ImmutableList;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;


@CommandLineProgramProperties(usage = "Count bases in dataflow with filters and transformers", usageShort = "count bases", programGroup = DataFlowProgramGroup.class)
public class CountBasesFiltered extends DataflowReadsPipeline {

    public CountBasesFiltered() {
        super(null, ImmutableList.of(ReadFilterLibrary.MAPPED), ImmutableList.of());
    }
}
