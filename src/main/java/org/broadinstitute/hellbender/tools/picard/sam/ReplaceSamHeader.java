package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Replace the SAMFileHeader in a SAM file with the given header. " +
                "Validation is minimal.  It is up to the user to ensure that all the elements referred to in the SAMRecords " +
                "are present in the new header.  Sort order of the two input files must be the same.",
        usageShort = "Replace the SAMFileHeader in a SAM file with the given header",
        programGroup = ReadProgramGroup.class
)
public class ReplaceSamHeader extends PicardCommandLineProgram {

    @Argument(doc = "SAM file from which SAMRecords will be read.", shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "SAM file from which SAMFileHeader will be read.")
    public File HEADER;

    @Argument(doc = "SAMFileHeader from HEADER file will be written to this file, followed by SAMRecords from INPUT file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(HEADER);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SAMFileHeader replacementHeader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(HEADER);

        if (BamFileIoUtils.isBamFile(INPUT)) {
            blockCopyReheader(replacementHeader);
        } else {
            standardReheader(replacementHeader);
        }

        return null;
    }

    private void standardReheader(final SAMFileHeader replacementHeader) {
        final SamReader recordReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).validationStringency(ValidationStringency.SILENT).open(INPUT);
        if (replacementHeader.getSortOrder() != recordReader.getFileHeader().getSortOrder()) {
            throw new UserException("Sort orders of INPUT (" + recordReader.getFileHeader().getSortOrder().name() +
                    ") and HEADER (" + replacementHeader.getSortOrder().name() + ") do not agree.");
        }
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(replacementHeader, true, OUTPUT);

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(ReplaceSamHeader.class));
        for (final SAMRecord rec : recordReader) {
            rec.setHeader(replacementHeader);
            writer.addAlignment(rec);
            progress.record(rec);
        }
        writer.close();
        CloserUtil.close(recordReader);
    }

    private void blockCopyReheader(final SAMFileHeader replacementHeader) {
        BamFileIoUtils.reheaderBamFile(replacementHeader, INPUT, OUTPUT, CREATE_MD5_FILE, CREATE_INDEX);
    }
}
