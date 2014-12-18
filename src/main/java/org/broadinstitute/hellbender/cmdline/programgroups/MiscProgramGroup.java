package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public class MiscProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Miscellaneous Tools"; }
    @Override
    public String getDescription() { return "A set of miscellaneous tools."; }
}
