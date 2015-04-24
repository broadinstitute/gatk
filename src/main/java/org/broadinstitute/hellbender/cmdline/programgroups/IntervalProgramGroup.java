package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

/**
 * Program group for tools that process intervals and associated overlapping records
 */
public final class IntervalProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() {
        return "Intervals";
    }

    @Override
    public String getDescription() {
        return "Tools for processing intervals and associated overlapping records";
    }

}
