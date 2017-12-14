package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that align reads, flag duplicates and recalibrate base qualities
 *
 * This class is a temporary PLACEHOLDER to be used only until the corresponding Picard group is available.
 */
public class BAMPreprocessingProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_BAM_PREPROCESSING; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_BAM_PREPROCESSING_SUMMARY; }
}