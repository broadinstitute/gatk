package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Program group for Example programs
 */
public class ExampleProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() {
        return HelpConstants.DOC_CAT_EXAMPLE;
    }

    @Override
    public String getDescription() {
        return HelpConstants.DOC_CAT_EXAMPLE_SUMMARY;
    }
}
