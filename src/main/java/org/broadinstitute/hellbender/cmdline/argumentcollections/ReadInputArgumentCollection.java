package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;
import java.util.List;


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

    @ArgumentCollection
    public ReadFilterArgumentCollection readFilterArgumentCollection = new ReadFilterArgumentCollection();

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line
     */
    public abstract List<File> getReadFiles();

    /**
     * Get the list of BAM/SAM/CRAM filenames specified at the command line
     */
    public abstract List<String> getReadFilesNames();

    /**
     * Get the read validation stringency specified at the command line, or the default value if none was specified
     * at the command line.
     */
    public ValidationStringency getReadValidationStringency() { return readValidationStringency; };

    public ReadFilterArgumentCollection getReadFilterArgumentCollection() {
        return readFilterArgumentCollection;
    }
}
