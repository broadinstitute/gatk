package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that count coverage, e.g. depth per allele
 */
public class CoverageAnalysisProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_COVERAGE_ANALYSIS; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_COVERAGE_ANALYSIS_SUMMARY; }
}
