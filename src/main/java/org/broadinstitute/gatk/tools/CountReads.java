package org.broadinstitute.gatk.tools;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
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
	usage = "Count reads.",
	usageShort = "Count reads",
        programGroup = SamOrBam.class
)
public class CountReads extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM file to count reads.")
    public File INPUT;

    public static void main(final String[] args) {
        new CountReads().instanceMain(args);
    }

    @Override
    protected int doWork() {
        long count = countReads();
        System.out.println(count + " reads");
        return 0;
    }

    /**
     * This is factored out of doWork only for unit testing.
     */
    long countReads() {
        long count=0;
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        for (final SAMRecord rec : in) {
            count++;
        }
        CloserUtil.close(in);
        return count;
    }
}