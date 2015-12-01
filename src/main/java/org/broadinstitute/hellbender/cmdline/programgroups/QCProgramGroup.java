package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public final class QCProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return "QC"; }
    @Override
    public String getDescription() { return "Tools for Diagnostics and Quality Control"; }
}
