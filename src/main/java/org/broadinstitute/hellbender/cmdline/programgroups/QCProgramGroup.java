package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public class QCProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return "QC"; }
    @Override
    public String getDescription() { return "Diagnostics and Quality Control Tools"; }
}
