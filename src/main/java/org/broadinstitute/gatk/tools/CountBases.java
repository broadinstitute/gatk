package org.broadinstitute.gatk.tools;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;

@CommandLineProgramProperties(
	usage = "Walks over the input data set, calculating the number of bases seen for diagnostic purposes.",
	usageShort = "Count bases",
        programGroup = SamOrBam.class
)
public class CountBases extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM file to count bases.")
    public File INPUT;

    public static void main(final String[] args) {
        new CountBases().instanceMain(args);
    }

    @Override
    protected int doWork() {
        long count = countBases();
        System.out.println(count + " bases");
        return 0;
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