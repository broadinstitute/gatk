package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.util.List;

/**
 * A tool to add comments to a BAM or CRAM file header. Effectively copies the file except for the addition of the @CO records
 * in the header. This tool does not support SAM files.
 *
 * @author jgentry
 */
@CommandLineProgramProperties(
        summary = "Adds one or more comments to the header of a specified BAM or CRAM file. Copies the file with the " +
                "modified header to a specified output file. Note that a block copying method is used to ensure efficient transfer to the " +
                "output file. SAM files are not supported",
        oneLineSummary = "Adds comments to the header of a BAM file",
        programGroup = ReadProgramGroup.class
)
public final class AddCommentsToBam extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input BAM/CRAM file to add a comment to the header")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output BAM file to write results")
    public File OUTPUT;

    @Argument(shortName = "C", doc = "Comments to add to the BAM file")
    public List<String> COMMENT;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        if (INPUT.getAbsolutePath().endsWith(".sam")) {
            throw new UserException("SAM files are not supported");
        }

        final SAMFileHeader samFileHeader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(INPUT);
        for (final String comment : COMMENT) {
            if (comment.contains("\n")) {
                throw new UserException("Comments can not contain a new line");
            }
            samFileHeader.addComment(comment);
        }

        BamFileIoUtils.reheaderBamFile(samFileHeader, INPUT, OUTPUT, CREATE_MD5_FILE, CREATE_INDEX);

        return null;
    }
}
