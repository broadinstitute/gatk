package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that analyze read coverage to detect copy number variants
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CopyNumberProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_CNV; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_CNV_SUMMARY; }
}
