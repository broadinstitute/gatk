package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.utils.read.SamComparison;

import java.io.File;
import java.util.*;

/**
 * Wrapper CLP for SamComparison.
 */
@CommandLineProgramProperties(
        usage = "USAGE: CompareSAMs <SAMFile1> <SAMFile2>\n" +
                "Compares the headers of the two input SAM or BAM files, and, if possible, the SAMRecords. " +
                "For SAMRecords, compares only the readUnmapped flag, reference name, start position and strand. " +
                "Reports the number of SAMRecords that match, differ in alignment, are mapped in only one input, " +
                "or are missing in one of the files",
        usageShort = "Compares two input SAM or BAM files",
        programGroup = ReadProgramGroup.class
)
public class CompareSAMs extends PicardCommandLineProgram {

    @PositionalArguments(minElements = 2, maxElements = 2)
    public List<File> samFiles;

    @Override
    protected Object doWork() {
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        SamReader sam1 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(0));
        SamReader sam2 = factory.referenceSequence(REFERENCE_SEQUENCE).open(samFiles.get(1));
        SamComparison comparison = new SamComparison(sam1, sam2);
        comparison.printReport();
        if (comparison.areEqual()) {
            System.out.println("SAM files match.");
        } else {
            System.out.println("SAM files differ.");
        }
        CloserUtil.close(sam1);
        CloserUtil.close(sam2);
        return null;
    }
}
