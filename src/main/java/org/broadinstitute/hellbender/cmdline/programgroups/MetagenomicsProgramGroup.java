package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

/**
 * Tools that perform metagenomic analysis, e.g. microbial community composition and pathogen detection
 */
public class MetagenomicsProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_METAGENOMICS; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_METAGENOMICS_SUMMARY; }
}