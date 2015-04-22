package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public final class IlluminaProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Illumina"; }
    @Override
    public String getDescription() { return "Tools for parsing + manipulating raw data from Illumina sequencers."; }
}
