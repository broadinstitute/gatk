package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.TaggedInputFileArgument;

import java.util.List;

/**
 * An argument collection for use with tools that accept one or more input files containing reads
 * (eg., BAM/SAM/CRAM files), and require at least one such input.
 */
public final class RequiredReadInputArgumentCollection extends ReadInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
              doc = "BAM/SAM/CRAM file containing reads",
              optional = false)
    private List<TaggedInputFileArgument> readInputs;

    @Override
    public List<TaggedInputFileArgument> getReadInputs() { return readInputs; }
}
