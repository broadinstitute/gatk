package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;
import java.util.List;

/**
 * An argument collection for use with tools that accept zero or more input files containing reads
 * (eg., BAM/SAM/CRAM files).
 */
public class OptionalReadInputArgumentCollection implements ArgumentCollectionDefinition {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "One or more BAM/SAM/CRAM files containing reads", optional = true)
    public List<File> readFiles;

}
