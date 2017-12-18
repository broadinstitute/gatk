package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that manipulate variant call format (VCF) data
 *
 * This class is a temporary PLACEHOLDER to be used only until the corresponding Picard group is available.
 */
public class VCFManipulationProgramGroup  implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_VCF_MANIPULATION; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_VCF_MANIPULATION_SUMMARY; }
}
