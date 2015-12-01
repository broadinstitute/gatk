package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.read.mergealignment.AbstractAlignmentMerger;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = CleanSam.USAGE,
        oneLineSummary = CleanSam.USAGE,
        programGroup = ReadProgramGroup.class
)
public final class CleanSam extends PicardCommandLineProgram {

    static final String USAGE = "Cleans the provided SAM/BAM/CRAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads";

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input SAM/BAM/CRAM to be cleaned.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Where to write cleaned SAM/BAM/CRAM.")
    public File OUTPUT;

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     */
    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReaderFactory factory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);
        if (VALIDATION_STRINGENCY == ValidationStringency.STRICT) {
            factory.validationStringency(ValidationStringency.LENIENT);
        }
        final SamReader reader = factory.open(INPUT);
        try (final SAMFileWriter writer = createSAMWriter(OUTPUT, REFERENCE_SEQUENCE, reader.getFileHeader(), true);
             final CloseableIterator<SAMRecord> it = reader.iterator()) {

            final ProgressLogger progress = new ProgressLogger(logger);

            // If the read (or its mate) maps off the end of the alignment, clip it
            while (it.hasNext()) {
                final SAMRecord rec = it.next();

                // If the read (or its mate) maps off the end of the alignment, clip it
                AbstractAlignmentMerger.createNewCigarsIfMapsOffEndOfReference(rec);

                // check the read's mapping quality
                if (rec.getReadUnmappedFlag() && 0 != rec.getMappingQuality()) {
                    rec.setMappingQuality(0);
                }

                writer.addAlignment(rec);
                progress.record(rec);
            }
        }
        finally {
            CloserUtil.close(reader);
        }
        return null;
    }
}
