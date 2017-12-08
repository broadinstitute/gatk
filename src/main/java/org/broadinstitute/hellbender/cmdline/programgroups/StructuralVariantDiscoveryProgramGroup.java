package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that detect structural variants
 */
public class StructuralVariantDiscoveryProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_SV_DISCOVERY; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_SV_DISCOVERY_SUMMARY; }
}
