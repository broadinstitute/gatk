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
import java.text.DecimalFormat;
import java.text.NumberFormat;

@CommandLineProgramProperties(
	usage = "Walks over all input data, accumulating statistics such as total number of read\n" +
            "reads with QC failure flag set, number of duplicates, percentage mapped, etc.",
	usageShort = "A reimplementation of the 'samtools flagstat' subcommand in the GATK.",
        programGroup = SamOrBam.class
)
public class FlagStat extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM file.")
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
            StringBuilder builder = new StringBuilder(ret);
            NumberFormat percentFormatter = new DecimalFormat("#0.00");
            builder.append(readCount);
            builder.append(" in total\n");

            builder.append(QC_failure);
            builder.append(" QC failure\n");


            builder.append(duplicates);
            builder.append(" duplicates\n");

            builder.append(mapped);
            builder.append(" mapped (");
            builder.append(percentFormatter.format(( (float)mapped / (float)readCount ) * 100.0));
            builder.append("%)\n");

            builder.append(paired_in_sequencing);
            builder.append(" paired in sequencing\n");


            builder.append(read1);
            builder.append(" read1\n");

            builder.append(read2);
            builder.append(" read2\n");

            builder.append(properly_paired);
            builder.append(" properly paired (");
            builder.append(percentFormatter.format(( (float)properly_paired / (float)readCount ) * 100.0));
            builder.append("%)\n");


            builder.append(with_itself_and_mate_mapped);
            builder.append(" with itself and mate mapped\n");

            builder.append(singletons);
            builder.append(" singletons (");
            builder.append(percentFormatter.format(( (float)singletons / (float)readCount ) * 100.0));
            builder.append("%)\n");


            builder.append(with_mate_mapped_to_a_different_chr);
            builder.append(" with mate mapped to a different chr\n");

            builder.append(with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5);
            builder.append(" with mate mapped to a different chr (mapQ>=5)");

            return builder.toString();
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

