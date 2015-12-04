package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

/**
 * Copy Number tool group.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CopyNumberProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "Copy Number Analysis";
    }

    @Override
    public String getDescription() {
        return "Tools to analyze copy number data.";
    }
}
