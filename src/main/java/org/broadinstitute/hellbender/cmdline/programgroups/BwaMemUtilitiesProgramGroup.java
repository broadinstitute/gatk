package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Program group for using misc part of JNI binding to Bwa Mem
 *
 * This Program Group is OBSOLETE and references to it should be replaced with a reference to an approved group.
 */
public final class BwaMemUtilitiesProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_BWA_MEM_UTILS; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_BWA_MEM_UTILS_SUMMARY; }
}
