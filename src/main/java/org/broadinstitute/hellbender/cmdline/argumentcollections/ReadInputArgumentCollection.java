package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.Serializable;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;


/**
 * An abstract argument collection for use with tools that accept input files containing reads
 * (eg., BAM/SAM/CRAM files).
 */
public abstract class ReadInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME,
            doc = "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true,
            optional=true)
    protected ValidationStringency readValidationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Argument(fullName = StandardArgumentDefinitions.READ_INDEX_LONG_NAME, shortName = StandardArgumentDefinitions.READ_INDEX_SHORT_NAME,
            doc = "Indices to use for the read inputs. If specified, an index must be provided for every read input " +
                    "and in the same order as the read inputs. If this argument is not specified, the path to the index " +
                    "for each input will be inferred automatically.",
            common = true,
            optional = true)
    protected List<GATKPath> readIndices;

    /**
     * Get the list of BAM/SAM/CRAM inputs specified at the command line.
     * GATKPath is the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */
    public abstract List<GATKPath> getReadPathSpecifiers();

    /**
     * Get the list of BAM/SAM/CRAM inputs specified at the command line.
     * GATKPath is the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */

    public List<Path> getReadPaths() {
        return getReadPathSpecifiers().stream().map(GATKPath::toPath).collect(Collectors.toList());
    }

    /**
     * @return The list of indices to be used with the read inputs, or {@code null} if none were specified and the indices should be
     *         inferred automatically.
     *
     *         If explicit indices are specified, they must be specified for all read inputs, and are assumed to be in the same
     *         order as the read inputs.
     */
    public List<Path> getReadIndexPaths() {
        if ( readIndices == null || readIndices.isEmpty() ) {
            return null;
        }

        return readIndices.stream().map(GATKPath::toPath).collect(Collectors.toList());
    }

    /**
     * Get the read validation stringency specified at the command line, or the default value if none was specified
     * at the command line.
     */
    public ValidationStringency getReadValidationStringency() { return readValidationStringency; };
}
