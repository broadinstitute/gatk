package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;

/**
 * Converts a BAM file to human-readable SAM output or vice versa
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = "Convert a SAM/BAM/CRAM file to a SAM/BAM/CRAM file (i.e., changes the format).\n" +
                "Input and output formats are determined by file extension.",
        oneLineSummary = "Convert a SAM/BAM/CRAM file to a SAM/BAM/CRAM file",
        programGroup = ReadProgramGroup.class
)
public final class SamFormatConverter extends PicardCommandLineProgram {

    @Argument(doc = "The input SAM/BAM/CRAM file.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "The ouput SAM/BAM/CRAM file. ",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        try (final SAMFileWriter writer = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, reader.getFileHeader(), true)) {

            final ProgressLogger progress = new ProgressLogger(logger);
            for (final SAMRecord rec : reader) {
                writer.addAlignment(rec);
                progress.record(rec);
            }
        }
        CloserUtil.close(reader);
        return null;
    }
}
