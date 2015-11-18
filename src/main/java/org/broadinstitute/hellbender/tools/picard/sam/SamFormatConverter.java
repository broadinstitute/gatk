package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;

/**
 * Converts a BAM file to human-readable SAM output or vice versa
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = "Convert a BAM file to a SAM file, or SAM to BAM.\n" +
                "Input and output formats are determined by file extension.",
        oneLineSummary = "Convert a BAM file to a SAM file, or a SAM to a BAM",
        programGroup = ReadProgramGroup.class
)
public final class SamFormatConverter extends PicardCommandLineProgram {
    @Argument(doc = "The BAM or SAM file to parse.", shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    @Argument(doc = "The BAM or SAM output file. ", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
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
