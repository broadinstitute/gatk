package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

public final class FastaProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_FASTA; }
    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_FASTA_SUMMARY; }
}
