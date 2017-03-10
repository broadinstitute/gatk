package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Program group for using non-essential part of JNI binding to Bwa Mem
 */
public final class BwaMemUtilitiesProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_BWA_MEM_UTILS; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_BWA_MEM_UTILS_SUMMARY; }
}
