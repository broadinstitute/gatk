package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public final class SparkPipelineProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Spark pipelines";
    }

    @Override
    public String getDescription() {
        return "Pipelines that combine tools and use Apache Spark for scaling out (experimental)";
    }
}

