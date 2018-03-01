package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

/**
 * Documentation-only arguments for GATK.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @implNote this class is not final to allow downstream projects to extend documentation-only arguments, but keep GATK
 * ones.
 */
public class GATKDocOnlyArgumentCollection implements DocOnlyArgumentCollection {

    // This option is here for documentation completeness.
    // This is actually parsed out in Main to initialize configuration files because
    // we need to have the configuration completely set up before we create our CommandLinePrograms.
    // (Some of the CommandLinePrograms have default values set to config values, and these are loaded
    // at class load time as static initializers).
    @Argument(fullName = StandardArgumentDefinitions.GATK_CONFIG_FILE_OPTION,
            doc = "A configuration file to use with the GATK.",
            common = true,
            optional = true)
    public String GATK_CONFIG_FILE = null;

}
