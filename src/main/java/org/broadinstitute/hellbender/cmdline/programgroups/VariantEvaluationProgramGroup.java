package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that evaluate and refine variant calls, e.g. by adding annotations that are not calculated during variant calling
 *
 * This class is a temporary PLACEHOLDER to be used only until the corresponding Picard group is available.
 */
public class VariantEvaluationProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_VARIANT_EVALUATION; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_VARIANT_EVALUATION_SUMMARY; }
}