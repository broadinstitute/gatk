package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.Option;
import org.broadinstitute.hellbender.cmdline.StandardOptionDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;

import java.io.File;

@CommandLineProgramProperties(
	usage = "Walks over the input data set, calculating the number of bases seen for diagnostic purposes.",
	usageShort = "Count bases",
        programGroup = ReadProgramGroup.class
)
public class CountBases extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM file to count bases.")
    public File INPUT;

    @Override
    protected Long doWork() {
        final long count = countBases();
        System.out.println(count + " bases");
        return count;
    }

    /**
     * This is factored out of doWork only for unit testing.
     */
    long countBases() {
        long count=0;
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        for (final SAMRecord rec : in) {
            count += rec.getReadLength();
        }
        CloserUtil.close(in);
        return count;
    }
}