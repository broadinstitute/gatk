package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.PositionalArguments;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.read.SamComparison;

import java.io.File;
import java.util.List;

/**
 * Wrapper CLP for SamComparison.
 */
@CommandLineProgramProperties(
        summary = "USAGE: CompareSAMs <SAMFile1> <SAMFile2>\n" +
                "Compares the headers of the two input SAM/BAM/CRAM files, and, if possible, the SAMRecords. " +
                "For SAMRecords, compares only the readUnmapped flag, reference name, start position and strand. " +
                "Reports the number of SAMRecords that match, differ in alignment, are mapped in only one input, " +
                "or are missing in one of the files",
        oneLineSummary = "Compares two input SAM/BAM/CRAM files",
        programGroup = ReadProgramGroup.class
)
public final class CompareSAMs extends PicardCommandLineProgram {

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> samFiles;

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(samFiles.get(0));
        IOUtil.assertFileIsReadable(samFiles.get(1));

        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(VALIDATION_STRINGENCY);
        SamReader sam1 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(0));
        SamReader sam2 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(1));
        SamComparison comparison = new SamComparison(sam1, sam2);
        comparison.printReport();
        if (comparison.areEqual()) {
            System.out.println("Files match.");
        } else {
            System.out.println("Files differ.");
        }
        CloserUtil.close(sam1);
        CloserUtil.close(sam2);
        return comparison.areEqual();
    }
}
