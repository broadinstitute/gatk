package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

/**
 * Program group for tools that manipulate variants and associated metadata
 */
public final class VariantProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "VCF";
    }

    @Override
    public String getDescription() {
        return "Tools for manipulating variants and associated metadata";
    }
}
