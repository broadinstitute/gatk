package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Program group for tools that process intervals and associated overlapping records
 */
public final class IntervalProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_INTERVALS; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_INTERVALS_SUMMARY; }

}
