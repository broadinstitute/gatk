package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.hadoop.fs.FileSystem;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.InputStream;
import java.net.URI;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.WeakHashMap;

/**
 * Represents a collection of assembly files  of file in the assembly collection
 */
public class AssemblyCollection {

    private final Path fastqDir;
    private final String fastqFileFormat;

    private final WeakHashMap<String, List<Template>> inputsById;

    public AssemblyCollection(final String fastqDir, final String fastqFileFormat) {
        Utils.nonNull(fastqDir);
        this.fastqDir = FileSystems.getFileSystem(URI.create(fastqDir)).getPath(fastqDir);
        this.fastqFileFormat = Utils.nonNull(fastqFileFormat);
        this.inputsById = new WeakHashMap<>(10);
    }

    public List<Template> inputs(final int assemblyNumber) {
        Utils.nonNull(assemblyNumber);
        final Path path = fastqDir.resolve(String.format(fastqFileFormat, assemblyNumber));
        if (!Files.exists(path)) {
            throw new UserException.CouldNotReadInput(path, "missing input file for assembly number " + assemblyNumber);
        } else {
            final List<SVFastqUtils.FastqRead> fastqReads = SVFastqUtils.readFastqFile(path.toString());
            final int numberOfFastqReads = fastqReads.size();
            if ((numberOfFastqReads & 1) == 1) {
                throw new UserException.BadInput("There is a odd number of reads in the assembly fastq file: " + path);
            }
            for (int i = 0; i < fastqReads.size(); i += 2) {
                final SVFastqUtils.FastqRead first = fastqReads.get(i);
                final SVFastqUtils.FastqRead second = fastqReads.get(i + 1);
                if (!first.getName().equals(second.getName())) {
                    throw new UserException.BadInput(
                            String.format("consecutive reads don't have the same name '%s' != '%s' in '%s'", first.getName(), second.getName(), path));
                } else if (first.getFragmentNumber().orElseThrow(() -> new UserException.BadInput("missing first fragment for")) != 1) {
                    throw new UserException.BadInput(
                            String.format("first in pair fragment number is not 1 (%d) for '%s' in '%s'", first.getFragmentNumber().getAsInt(), first.getName(), path));
                } else if (second.getFragmentNumber().orElseThrow(() -> new UserException.BadInput("missing first fragment for")) != 2) {
                    throw new UserException.BadInput(
                            String.format("first in pair fragment number is not 2 (%d) for '%s' in '%s'", second.getFragmentNumber().getAsInt(), first.getName(), path));
                } else {
                    Template.create(first.getName(), Arrays.asList(first, second, r -> new Template.Fragment(r.getName(), r.getFragmentNumber)))
                }
            }
            Template.create()
        }
    }

}
