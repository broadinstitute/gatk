package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = "Sorts the input SAM or BAM.\n" +
                "Input and output formats are determined by file extension.",
        oneLineSummary = "Sorts a SAM or BAM file",
        programGroup = ReadProgramGroup.class
)
public final class SortSam extends PicardCommandLineProgram {

    @Argument(doc = "The BAM or SAM file to sort.", shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "The sorted BAM or SAM output file. ", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(shortName = StandardArgumentDefinitions.SORT_ORDER_SHORT_NAME, doc = "Sort order of output file")
    public SAMFileHeader.SortOrder SORT_ORDER;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        ;
        reader.getFileHeader().setSortOrder(SORT_ORDER);
        try (final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, OUTPUT)) {
            writer.setProgressLogger(
                    new ProgressLogger(logger, (int) 1e7, "Wrote", "records from a sorting collection"));

            final ProgressLogger progress = new ProgressLogger(logger, (int) 1e7, "Read");
            for (final SAMRecord rec : reader) {
                writer.addAlignment(rec);
                progress.record(rec);
            }

            logger.info("Finished reading inputs, merging and writing to output now.");

        }
        CloserUtil.close(reader);
        return null;
    }
}
