package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * An argument collection for use with tools that accept zero or more input files containing reads
 * (eg., BAM/SAM/CRAM files).
 */
public final class OptionalReadInputArgumentCollection extends ReadInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "BAM/SAM/CRAM file containing reads",
            optional = true,
            common = true)
    private List<GATKPath> readInputNames = new ArrayList<>();

    /**
     * Get the list of BAM/SAM/CRAM inputs specified at the command line.
     * GATKPath is the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */
    @Override
    public List<GATKPath> getRawReadPathSpecifiers() { return Collections.unmodifiableList(readInputNames); }
}
