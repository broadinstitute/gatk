package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public class DataFlowProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Dataflow tests";
    }

    @Override
    public String getDescription() {
        return "Programs to test out google data flow";
    }
}

