package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that process genomic intervals in various formats
 *
 * This class is a temporary PLACEHOLDER to be used only until the corresponding Picard group is available.
 */
public final class IntervalsProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_INTERVALS; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_INTERVALS_SUMMARY; }

}
