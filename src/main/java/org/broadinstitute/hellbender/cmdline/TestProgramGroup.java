package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * only for testing
 */
public final class TestProgramGroup implements CommandLineProgramGroup {
    @Override
    public String getName() { return HelpConstants.DOC_CAT_TEST; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_TEST_SUMMARY; }
}