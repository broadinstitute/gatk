package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

public final class FastaProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Fasta"; }
    @Override
    public String getDescription() { return "Tools for analysis and manipulation of files in fasta format"; }
}
