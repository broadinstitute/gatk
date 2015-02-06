package org.broadinstitute.hellbender.tools.picard;

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
 * Converts a BAM file to human-readable SAM output or vice versa
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Convert a BAM file to a SAM file, or SAM to BAM.\n" +
                "Input and output formats are determined by file extension.",
        usageShort = "Convert a BAM file to a SAM file, or a SAM to a BAM",
        programGroup = ReadProgramGroup.class
)
public class SamFormatConverter extends PicardCommandLineProgram {
    @Argument(doc = "The BAM or SAM file to parse.", shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    @Argument(doc = "The BAM or SAM output file. ", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(reader.getFileHeader(), true, OUTPUT, REFERENCE_SEQUENCE);

        if (CREATE_INDEX && writer.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new UserException("Can't CREATE_INDEX unless sort order is coordinate");
        }

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(SamFormatConverter.class));
        for (final SAMRecord rec : reader) {
            writer.addAlignment(rec);
            progress.record(rec);
        }
        CloserUtil.close(reader);
        writer.close();
        return null;
    }
}
