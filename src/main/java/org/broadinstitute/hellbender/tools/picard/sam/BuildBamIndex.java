package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Command line program to generate a BAM index (.bai) file from a BAM (.bam) file
 *
 * @author Martha Borkan
 */
@CommandLineProgramProperties(
        summary = "Generates a BAM index (.bai) file.",
        oneLineSummary = "Generates a BAM index (.bai) file",
        programGroup = ReadProgramGroup.class
)
public final class BuildBamIndex extends PicardCommandLineProgram {

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "A BAM file or URL to process. Must be sorted in coordinate order.")
    public String INPUT;

    URL inputUrl = null;   // INPUT as URL
    File inputFile = null; // INPUT as File, if it can't be interpreted as a valid URL

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.\n" +
                    "If input is a URL and OUTPUT is unspecified, defaults to a file in the current directory.", optional = true)
    public File OUTPUT;

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then iterates through
     * all the records generating a BAM Index, then writes the bai file.
     */
    protected Object doWork() {
        try {
            inputUrl = new URL(INPUT);
        } catch (MalformedURLException e) {
            inputFile = new File(INPUT);
        }

        // set default output file - input-file.bai
        if (OUTPUT == null) {

            final String baseFileName;
            if (inputUrl != null) {
                final String path = inputUrl.getPath();
                final int lastSlash = path.lastIndexOf("/");
                baseFileName = path.substring(lastSlash + 1, path.length());
            } else {
                baseFileName = inputFile.getAbsolutePath();
            }

            if (baseFileName.endsWith(BamFileIoUtils.BAM_FILE_EXTENSION)) {

                final int index = baseFileName.lastIndexOf(".");
                OUTPUT = new File(baseFileName.substring(0, index) + BAMIndex.BAMIndexSuffix);

            } else {
                OUTPUT = new File(baseFileName + BAMIndex.BAMIndexSuffix);
            }
        }

        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader bam;

        if (inputUrl != null) {
            // remote input
            bam = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE)
                    .disable(SamReaderFactory.Option.EAGERLY_DECODE)
                    .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                    .open(SamInputResource.of(inputUrl));
        } else {
            // input from a normal file
            IOUtil.assertFileIsReadable(inputFile);
            bam = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE)
                    .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                    .open(inputFile);
        }

        if (bam.type() != SamReader.Type.BAM_TYPE) {
            throw new SAMException("Input file must be bam file, not sam file.");
        }

        if (!bam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
            throw new SAMException("Input bam file must be sorted by coordinates");
        }

        BAMIndexer.createIndex(bam, OUTPUT);

        logger.info("Successfully wrote bam index file " + OUTPUT);
        CloserUtil.close(bam);
        return null;
    }
}
