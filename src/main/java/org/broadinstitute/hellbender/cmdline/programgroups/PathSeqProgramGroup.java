package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools for detecting pathogens using the PathSeq pipeline.
 */
public final class PathSeqProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_SPARK_PATHSEQ; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_SPARK_PATHSEQ_SUMMARY; }
}
