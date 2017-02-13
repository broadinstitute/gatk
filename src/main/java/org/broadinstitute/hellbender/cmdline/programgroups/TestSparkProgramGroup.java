package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

public class TestSparkProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_SPARK_VALIDATION; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_SPARK_VALIDATION_SUMMARY; }
}
