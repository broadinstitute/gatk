package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.File;
import java.io.Serializable;

/**
 * Collection of arguments that control how we create output alignments.
 */
public class OutputAlignmentArgumentCollection implements Serializable {

    private static final long serialVersionUID = 1L;

    private final ReferenceInputArgumentCollection referenceArguments;

    /**
     * Creates a new output alignment argument collection given the reference arguments.
     * @param referenceArguments the reference arguments, must not be {@code null}.
     */
    public OutputAlignmentArgumentCollection(final ReferenceInputArgumentCollection referenceArguments) {
        this.referenceArguments = Utils.nonNull(referenceArguments);
    }

    @Argument(fullName= StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_SHORT_NAME,
            doc = "If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.", optional=true, common = true)
    public boolean createOutputBamIndex = true;

    @Argument(fullName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_LONG_NAME,
            shortName=StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_SHORT_NAME,
            doc = "If true, create a MD5 digest for any BAM/SAM/CRAM file created", optional=true, common = true)
    public boolean createOutputBamMD5 = false;

    /*
     * Create a common SAMFileWriter using the reference and read header for this tool.
     *
     * @param outputFile    - if this file has a .cram extension then a reference is required. Can not be null.
     * @param preSorted     - if true then the records must already be sorted to match the header sort order
     *
     * @throws UserException if outputFile ends with ".cram" and no reference is provided
     * @return SAMFileWriter
     */
    public final SAMFileGATKReadWriter createSAMWriter(final File outputFile, final boolean preSorted, final SAMFileHeader header) {
        if (!referenceArguments.hasReference() && IOUtils.isCramFile(outputFile)) {
            throw new UserException.MissingReference("A reference file is required for writing CRAM files");
        }

        return new SAMFileGATKReadWriter(
                ReadUtils.createCommonSAMWriter(
                        outputFile,
                        referenceArguments.getReferenceFile(),
                        header,
                        preSorted,
                        createOutputBamIndex,
                        createOutputBamMD5
                )
        );
    }

}
