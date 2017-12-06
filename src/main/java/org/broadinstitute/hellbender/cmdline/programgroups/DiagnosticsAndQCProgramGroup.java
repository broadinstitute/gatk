package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that collect sequencing quality-related and comparative metrics
 *
 * This class is a temporary PLACEHOLDER to be used only until the corresponding Picard group is available.
 */
public class DiagnosticsAndQCProgramGroup  implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_DIAGNOSTICS_AND_QC; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_DIAGNOSTICS_AND_QC_SUMMARY; }
}