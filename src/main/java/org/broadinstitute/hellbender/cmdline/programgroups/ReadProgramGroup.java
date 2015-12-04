package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public final class ReadProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return "SAM/BAM/CRAM"; }
    @Override
    public String getDescription() { return "Tools for manipulating read-level data (SAM/BAM/CRAM)"; }
}
