package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public final class StructuralVariationProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Spark tools for structural variation analysis";
    }

    @Override
    public String getDescription() {
        return "Structural variation analysis tools that use Apache Spark for scaling out (experimental)";
    }
}
