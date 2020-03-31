package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An argument collection for use with tools that accept one or more input files containing reads
 * (eg., BAM/SAM/CRAM files), and require at least one such input.
 */
public final class RequiredReadInputArgumentCollection extends ReadInputArgumentCollection {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "BAM/SAM/CRAM file containing reads", optional = false, common = true)
    public List<GATKPathSpecifier> readFilesNames;

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line.
     * Paths are the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */
    public List<GATKPathSpecifier> getReadPathSpecifiers() { return Collections.unmodifiableList(readFilesNames);}

    public List<Path> getReadPaths() {
        return getReadPathSpecifiers().stream().map(GATKPathSpecifier::toPath).collect(Collectors.toList());
    };
}
