package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = "Replace the SAMFileHeader in a SAM/BAM file with the given header. " +
                "Validation is minimal.  It is up to the user to ensure that all the elements referred to in the SAMRecords " +
                "are present in the new header.  Sort order of the two input files must be the same.",
        oneLineSummary = "Replace the SAMFileHeader in a SAM/BAM file with the given header",
        programGroup = ReadProgramGroup.class
)
public final class ReplaceSamHeader extends PicardCommandLineProgram {

    @Argument(doc = "SAM file from which SAMRecords will be read.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "SAM file from which SAMFileHeader will be read.")
    public File HEADER;

    @Argument(doc = "SAMFileHeader from HEADER file will be written to this file, followed by SAMRecords from INPUT file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
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
        try (final SAMFileWriter writer = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, replacementHeader, true)) {

            final ProgressLogger progress = new ProgressLogger(logger);
            for (final SAMRecord rec : recordReader) {
                rec.setHeaderStrict(replacementHeader);
                writer.addAlignment(rec);
                progress.record(rec);
            }
        }
        CloserUtil.close(recordReader);
    }

    private void blockCopyReheader(final SAMFileHeader replacementHeader) {
        BamFileIoUtils.reheaderBamFile(replacementHeader, INPUT, OUTPUT, CREATE_MD5_FILE, CREATE_INDEX);
    }
}
