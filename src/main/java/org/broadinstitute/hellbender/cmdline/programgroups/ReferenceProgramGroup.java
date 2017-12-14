package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that analyze and manipulate FASTA format references
 *
 * This class is a temporary PLACEHOLDER to be used only until the corresponding Picard group is available.
 */
public class ReferenceProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_REFERENCE; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_REFERENCE_SUMMARY; }
}
