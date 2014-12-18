package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public class Testing implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Unit Testing"; }
    @Override
    public String getDescription() { return "Unit testing"; }
}
