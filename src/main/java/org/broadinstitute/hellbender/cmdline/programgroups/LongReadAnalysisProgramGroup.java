package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

public class LongReadAnalysisProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_LR_ANALYSIS; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_LR_ANALYSIS_SUMMARY; }
}
