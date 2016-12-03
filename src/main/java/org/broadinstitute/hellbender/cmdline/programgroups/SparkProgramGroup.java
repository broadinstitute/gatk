package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public final class SparkProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Spark tools";
    }

    @Override
    public String getDescription() {
        return "Tools that use Apache Spark for scaling out (experimental)";
    }
}

