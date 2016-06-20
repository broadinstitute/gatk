package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.TaggedInputFileArgument;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;


/**
 * An abstract argument collection for use with tools that accept input files containing reads
 * (eg., BAM/SAM/CRAM files).
 */
public abstract class ReadInputArgumentCollection implements ArgumentCollectionDefinition {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME,
            doc = "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true,
            optional=true)
    public ValidationStringency readValidationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line
     */
    public List<File> getReadFiles() {
        return getReadInputs().stream().map(TaggedInputFileArgument::getFile).collect(Collectors.toList());
    }

    /**
     * Get the list of BAM/SAM/CRAM filenames specified at the command line
     */
    public List<String> getReadFileNames() {
        return getReadInputs().stream().map(TaggedInputFileArgument::getFilePath).collect(Collectors.toList());
    }

    /**
     * Returns a map of symbolicName -> ReadInput.
     * The map is sorted by symbolic name
     */
    public SortedMap<String, TaggedInputFileArgument> getInputsBySymbolicName(){
        final List<TaggedInputFileArgument> readInputs = getReadInputs();
        final SortedMap<String, TaggedInputFileArgument> result = new TreeMap<>();//sort by name
        readInputs.forEach(ri -> result.put(ri.getName(), ri));
        return result;
    }

    /**
     * Returns all ReadInputs in this argument collection.
     */
    public abstract List<TaggedInputFileArgument> getReadInputs();

    /**
     * Get the read validation stringency specified at the command line, or the default value if none was specified
     * at the command line.
     */
    public ValidationStringency getReadValidationStringency() { return readValidationStringency; }
}
