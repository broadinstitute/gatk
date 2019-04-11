package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that performs methylation calling and methylation-based coverage for bisulfite BAMs
 */
public class MethylationProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_METHYLATION_DISCOVERY; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_METHYLATION_DISCOVERY_SUMMARY; }
}
