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
import java.text.DecimalFormat;
import java.text.NumberFormat;

@CommandLineProgramProperties(
	usage = "Walks over all input data, accumulating statistics such as total number of read\n" +
            "reads with QC failure flag set, number of duplicates, percentage mapped, etc.",
	usageShort = "A reimplementation of the 'samtools flagstat' subcommand.",
        programGroup = ReadProgramGroup.class
)
public class FlagStat extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM or CRAM file.")
    public File INPUT;

    public static void main(final String[] args) {
        new FlagStat().instanceMain(args);
    }

    @Override
    protected int doWork() {
        FlagStatus sum = countReads();
        System.out.println(sum);
        return 0;
    }

    /**
     * This is factored out of doWork only for unit testing.
     */
    FlagStatus countReads() {
        FlagStatus sum = new FlagStatus();
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        for (final SAMRecord rec : in) {
            sum.add(rec);

        }
        CloserUtil.close(in);
        return sum;
    }

    // what comes out of the flagstat
    public final static class FlagStatus {
        long readCount = 0L;
        long QC_failure = 0L;
        long duplicates = 0L;
        long mapped = 0L;
        long paired_in_sequencing = 0L;
        long read1 = 0L;
        long read2 = 0L;
        long properly_paired = 0L;
        long with_itself_and_mate_mapped = 0L;
        long singletons = 0L;
        long with_mate_mapped_to_a_different_chr = 0L;
        long with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 = 0L;

        public String toString() {
            String ret = "";
            NumberFormat percentFormatter = new DecimalFormat("#0.00");

            return ret + readCount + " in total\n"
                       + QC_failure + " QC failure\n"
                       + duplicates + " duplicates\n"
                       + mapped + " mapped ("
                       + percentFormatter.format(((float) mapped / (float) readCount) * 100.0) + "%)\n"
                       + paired_in_sequencing + " paired in sequencing\n"
                       + read1 + " read1\n"
                       + read2 + " read2\n"
                       + properly_paired + " properly paired ("
                       + percentFormatter.format(((float) properly_paired / (float) readCount) * 100.0) + "%)\n"
                       + with_itself_and_mate_mapped + " with itself and mate mapped\n"
                       + singletons + " singletons ("
                       + percentFormatter.format(((float) singletons / (float) readCount) * 100.0) + "%)\n"
                       + with_mate_mapped_to_a_different_chr + " with mate mapped to a different chr\n"
                       + with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 + " with mate mapped to a different chr (mapQ>=5)";
        }

        public FlagStatus add(final SAMRecord read) {
            this.readCount++;

            if (read.getReadFailsVendorQualityCheckFlag()) {
                this.QC_failure++;
            }
            if (read.getDuplicateReadFlag()) {
                this.duplicates++;
            }
            if (!read.getReadUnmappedFlag()) {
                this.mapped++;
            }
            if (read.getReadPairedFlag()) {
                this.paired_in_sequencing++;

                if (read.getSecondOfPairFlag()) {
                    this.read2++;
                } else if (read.getReadPairedFlag()) {
                    this.read1++;
                }
                if (read.getProperPairFlag()) {
                    this.properly_paired++;
                }
                if (!read.getReadUnmappedFlag() && !read.getMateUnmappedFlag()) {
                    this.with_itself_and_mate_mapped++;

                    if (!read.getReferenceIndex().equals(read.getMateReferenceIndex())) {
                        this.with_mate_mapped_to_a_different_chr++;

                        if (read.getMappingQuality() >= 5) {
                            this.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5++;
                        }
                    }
                }
                if (!read.getReadUnmappedFlag() && read.getMateUnmappedFlag()) {
                    this.singletons++;
                }
            }

            return this;
        }
    }
}

