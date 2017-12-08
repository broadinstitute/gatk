package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that manipulate read data in SAM, BAM or CRAM format
 *
 * This class is a temporary PLACEHOLDER to be used only until the corresponding Picard group is available.
 */
public class ReadDataProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_READ_DATA; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_READ_DATA_SUMMARY; }
}