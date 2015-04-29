package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;

import java.io.File;
import java.util.List;

/**
 * An argument collection for use with tools that accept zero or more input files containing reads
 * (eg., BAM/SAM/CRAM files).
 */
public final class OptionalReadInputArgumentCollection extends ReadInputArgumentCollection {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "BAM/SAM/CRAM file containing reads", optional = true)
    public List<File> readFiles;

    @Override
    public List<File> getReadFiles() {
        return readFiles;
    }
}
