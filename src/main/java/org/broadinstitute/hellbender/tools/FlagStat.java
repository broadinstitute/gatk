package org.broadinstitute.hellbender.tools;

import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;

@CommandLineProgramProperties(
	usage = "Walks over all input data, accumulating statistics such as total number of read\n" +
            "reads with QC failure flag set, number of duplicates, percentage mapped, etc.",
	usageShort = "A reimplementation of the 'samtools flagstat' subcommand.",
    programGroup = ReadProgramGroup.class
)
public final class FlagStat extends ReadWalker {


    private FlagStatus sum = new FlagStatus();

    @Override
    public void apply( SAMRecord read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        sum.add(read);
    }

    @Override
    public Object onTraversalDone() {
        return sum;
    }

    // what comes out of the flagstat
    public final static class FlagStatus implements Serializable{
        private static final long serialVersionUID = 1L;
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

        public void merge(final FlagStatus that){
            this.readCount += that.readCount;
            this.QC_failure += that.QC_failure;
            this.duplicates += that.duplicates;
            this.mapped += that.mapped;
            this.paired_in_sequencing += that.paired_in_sequencing;
            this.read1 += that.read1;
            this.read2 += that.read2;
            this.properly_paired += that.properly_paired;
            this.with_itself_and_mate_mapped += that.with_itself_and_mate_mapped;
            this.singletons += that.singletons;
            this.with_mate_mapped_to_a_different_chr += that.with_mate_mapped_to_a_different_chr;
            this.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 += that.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5;

        }


        public FlagStatus add(final Read read){
            if(read == null){
                return this;
            }

            this.readCount++;

            if (read.getFailedVendorQualityChecks() != null && read.getFailedVendorQualityChecks()) {
                this.QC_failure++;
            }
            if (read.getDuplicateFragment() != null && read.getDuplicateFragment()) {
                this.duplicates++;
            }
            if (read.getAlignment() != null) {
                this.mapped++;
            }
            if (read.getNumberReads() != null && read.getNumberReads() == 2) {
                this.paired_in_sequencing++;

                if (read.getReadNumber() != null && read.getReadNumber() == 1) {
                    this.read2++;
                } else if (read.getReadNumber() != null && read.getReadNumber() == 0) {
                    this.read1++;
                }
                if (read.getProperPlacement() != null && read.getProperPlacement()) {
                    this.properly_paired++;
                }
                if (read.getAlignment() != null && read.getNextMatePosition() != null){
                    this.with_itself_and_mate_mapped++;

                    if (!read.getAlignment().getPosition().getReferenceName().equals(read.getNextMatePosition().getReferenceName())) {
                        this.with_mate_mapped_to_a_different_chr++;

                        if (read.getAlignment().getMappingQuality() >= 5) {
                            this.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5++;
                        }
                    }
                }
                if (read.getAlignment() != null && read.getNextMatePosition() == null) {
                    this.singletons++;
                }
            }
            return this;
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

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            FlagStatus that = (FlagStatus) o;

            if (QC_failure != that.QC_failure) return false;
            if (duplicates != that.duplicates) return false;
            if (mapped != that.mapped) return false;
            if (paired_in_sequencing != that.paired_in_sequencing) return false;
            if (properly_paired != that.properly_paired) return false;
            if (read1 != that.read1) return false;
            if (read2 != that.read2) return false;
            if (readCount != that.readCount) return false;
            if (singletons != that.singletons) return false;
            if (with_itself_and_mate_mapped != that.with_itself_and_mate_mapped) return false;
            if (with_mate_mapped_to_a_different_chr != that.with_mate_mapped_to_a_different_chr) return false;
            if (with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 != that.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5)
                return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = (int) (readCount ^ (readCount >>> 32));
            result = 31 * result + (int) (QC_failure ^ (QC_failure >>> 32));
            result = 31 * result + (int) (duplicates ^ (duplicates >>> 32));
            result = 31 * result + (int) (mapped ^ (mapped >>> 32));
            result = 31 * result + (int) (paired_in_sequencing ^ (paired_in_sequencing >>> 32));
            result = 31 * result + (int) (read1 ^ (read1 >>> 32));
            result = 31 * result + (int) (read2 ^ (read2 >>> 32));
            result = 31 * result + (int) (properly_paired ^ (properly_paired >>> 32));
            result = 31 * result + (int) (with_itself_and_mate_mapped ^ (with_itself_and_mate_mapped >>> 32));
            result = 31 * result + (int) (singletons ^ (singletons >>> 32));
            result = 31 * result + (int) (with_mate_mapped_to_a_different_chr ^ (with_mate_mapped_to_a_different_chr >>> 32));
            result = 31 * result + (int) (with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 ^ (with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 >>> 32));
            return result;
        }
    }
}

