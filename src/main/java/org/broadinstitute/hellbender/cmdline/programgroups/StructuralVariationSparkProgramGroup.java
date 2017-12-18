package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools for structural variation analysis that runs on Spark
 *
 * This Program Group is OBSOLETE and references to it should be replaced with a reference to an approved group.
 */
public final class StructuralVariationSparkProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_SPARK_SV; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_SPARK_SV_SUMMARY; }
}
