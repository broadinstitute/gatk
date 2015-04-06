package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.read.mergealignment.AbstractAlignmentMerger;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = CleanSam.USAGE,
        usageShort = CleanSam.USAGE,
        programGroup = ReadProgramGroup.class
)
public class CleanSam extends PicardCommandLineProgram {

    static final String USAGE = "Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads";

    @Argument(shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, doc = "Input SAM to be cleaned.")
    public File INPUT;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write cleaned SAM.")
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
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
        final CloseableIterator<SAMRecord> it = reader.iterator();
        final ProgressLogger progress = new ProgressLogger(Log.getInstance(CleanSam.class));

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

        writer.close();
        it.close();
        CloserUtil.close(reader);
        return null;
    }
}
