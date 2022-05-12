package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that perform variant calling and genotyping for short variants (SNPs, SNVs and Indels)
 */
public class FlowBasedProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_SHORT_FLOW_BASED; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_SHORT_FLOW_BASED_SUMMARY; }
}