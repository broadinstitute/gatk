package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An argument collection for use with tools that accept one or more input files containing reads
 * (eg., BAM/SAM/CRAM files), and require at least one such input.
 */
public final class RequiredReadInputArgumentCollection extends ReadInputArgumentCollection {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "BAM/SAM/CRAM file containing reads", optional = false)
    public List<String> readFilesNames;

    @Argument(fullName = "readIndex", shortName = "readIndex", doc = "BAM/CRAM indices", optional = true)
    private List<String> indices;

    @Override
    public List<File> getReadFiles() {
        ArrayList<File> ret = new ArrayList<>();
        for (String fn : readFilesNames) {
            ret.add(new File(fn));
        }
        return ret;
    }

    @Override
    public List<Path> getReadPaths() {
        ArrayList<Path> ret = new ArrayList<>();
        for (String fn : readFilesNames) {
            ret.add(IOUtils.getPath(fn));
        }
        return ret;
    }

    @Override
    public List<String> getReadFilesNames() {
        return new ArrayList<>(readFilesNames);
    }

    @Override
    public List<Path> getReadIndexPaths() {
        if ( indices == null || indices.isEmpty() ) {
            return null;
        }

        return indices.stream().map(index -> IOUtils.getPath(index)).collect(Collectors.toList());
    }
}
