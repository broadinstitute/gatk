package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public final class SparkProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Spark tests";
    }

    @Override
    public String getDescription() {
        return "Programs to test out Apache Spark";
    }
}

