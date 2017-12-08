package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that perform variant calling and genotyping for short variants (SNPs, SNVs and Indels)
 */
public class ShortVariantDiscoveryProgramGroup  implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_SHORT_VARIANT_DISCOVERY; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_SHORT_VARIANT_DISCOVERY_SUMMARY; }
}