package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.hellbender.cmdline.CommandLineProgramGroup;

/**
 * Exome analysis tool group.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ExomeAnalysisProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return "Exome Analysis";
    }

    @Override
    public String getDescription() {
        return "Tools to analyze exome data.";
    }
}
