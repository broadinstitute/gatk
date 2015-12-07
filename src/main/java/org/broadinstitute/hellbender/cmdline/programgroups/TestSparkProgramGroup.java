package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public class TestSparkProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Spark Validation tools";
    }

    @Override
    public String getDescription() {
        return "Tools written in Spark to compare aspects of two different files";
    }
}
