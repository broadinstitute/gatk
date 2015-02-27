package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

@CommandLineProgramProperties(
        usage = "Creates a hash code based on identifying information in the RG (read group) " +
                "records in a SAM file's header. This hash code changes any time read groups are added or removed " +
                "comparing one file's hash code to another tells you if the read groups in the BAM files are different.",
        usageShort = "Creates a hash code based on the read groups (RG) in the SAM or BAM header.",
        programGroup = ReadProgramGroup.class
)
public class CalculateReadGroupChecksum extends PicardCommandLineProgram {

    private static final String OUTPUT_FILE_EXTENSION = ".read_group_md5";

    @Argument(shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file. ")
    public File INPUT;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "The file to which the hash code should be written.", optional = true)
    public File OUTPUT;

    /**
     * Creates a file name (not including the path) for an RG MD5 file based on the name of the input file.
     */
    public static String getOutputFileName(final File inputFile) {
        return inputFile.getName() + OUTPUT_FILE_EXTENSION;
    }

    @Override
    protected Object doWork() {
        final File output = (OUTPUT == null) ? new File(INPUT.getParentFile(), getOutputFileName(INPUT)) : OUTPUT;

        IOUtil.assertFileIsWritable(output);
        final String hashText = SAMUtils.calculateReadGroupRecordChecksum(INPUT, REFERENCE_SEQUENCE);

        try {
            final FileWriter outputWriter = new FileWriter(output);
            outputWriter.write(hashText);
            outputWriter.close();
        } catch (final IOException ioe) {
            throw new UserException("Could not write the computed hash (" + hashText + ") to the output file: " + ioe.getMessage(), ioe);
        }
        return null;
    }
}
