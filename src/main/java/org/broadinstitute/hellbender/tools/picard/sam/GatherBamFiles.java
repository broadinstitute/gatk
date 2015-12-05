package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;
import java.util.List;

/**
 * Program to perform a rapid "gather" operation on BAM files after a scatter operations where
 * the same process has been performed on different regions of a BAM file creating many smaller
 * BAM files that now need to be concatenated back together.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Concatenates one or more BAM or CRAM files together as efficiently as possible. Assumes that the " +
                "list of files provided as INPUT are in the order that they should be concatenated. For BAM files, simply concatenates the bodies " +
                "of the files, retaining the header from the first file, via copying of the gzip blocks directly for speed. " +
                "Also supports generation of an MD5 on the output and indexing of the output file. Does not " +
                "support SAM files.",
        oneLineSummary = "Concatenates one or more BAM files together as efficiently as possible",
        programGroup = ReadProgramGroup.class
)
public final class GatherBamFiles extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "One or more BAM or CRAM files or text files containing lists of BAM files one per line.")
    public List<File> INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file to write.")
    public File OUTPUT;

    private static final Logger log = LogManager.getLogger();

    @Override
    protected Object doWork() {
        final List<File> inputs = IOUtil.unrollFiles(INPUT, BamFileIoUtils.BAM_FILE_EXTENSION, ".sam", CramIO.CRAM_FILE_EXTENSION);
        for (final File f : inputs) IOUtil.assertFileIsReadable(f);
        IOUtil.assertFileIsWritable(OUTPUT);

        if (determineBlockCopyingStatus(inputs)) {
            BamFileIoUtils.gatherWithBlockCopying(inputs, OUTPUT, CREATE_INDEX, CREATE_MD5_FILE);
        } else {
            gatherNormally(inputs, OUTPUT, REFERENCE_SEQUENCE);
        }

        return null;
    }

    private boolean determineBlockCopyingStatus(final List<File> inputs) {
        boolean useBlockCopying = true;
        for (final File f : inputs) {
            if (!BamFileIoUtils.isBamFile(f)) {
                useBlockCopying = false;
            }
        }
        return useBlockCopying;
    }

    /**
     * Simple implementation of a gather operations that uses SAMFileReaders and Writers in order to concatenate
     * multiple BAM files.
     */
    private void gatherNormally(final List<File> inputs, final File output, final File referenceFasta) {
        final SAMFileHeader header;
        {
            header = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).getFileHeader(inputs.get(0));
        }

        try (final SAMFileWriter out = createSAMWriter(output, referenceFasta, header, false)) {
            for (final File f : inputs) {
                log.info("Gathering " + f.getAbsolutePath());
                final SamReader in = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(f);
                for (final SAMRecord rec : in) {
                    out.addAlignment(rec);
                }
                CloserUtil.close(in);
            }
        }
    }

}
