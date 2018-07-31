package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Program group for use with internal test CommandLinePrograms only.
 */
public final class TestProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_TEST; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_TEST_SUMMARY; }
}