package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

/**
 * Tools for structural variation analysis that runs on Spark.
 */
public final class StructuralVariationSparkProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Spark tools for structural variation analysis";
    }

    @Override
    public String getDescription() {
        return "Structural variation analysis tools that use Apache Spark for scaling out (experimental)";
    }
}
