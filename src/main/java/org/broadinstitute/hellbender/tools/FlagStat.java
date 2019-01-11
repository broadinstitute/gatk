package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Accumulate flag statistics given a BAM file, e.g. total number of reads with QC failure flag set, number of
 * duplicates, percentage mapped etc.
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>A BAM file containing aligned read data</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>Accumulated flag statistics</li>
 * </ul>
 *
 * <h3>Example Usage</h3>
 * <pre>
 *   gatk FlagStat \
 *     -I input.bam
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
	summary = "Accumulate flag statistics given a BAM file, e.g. total number of reads with QC failure flag set, " +
            "number of duplicates, percentage mapped etc.",
	oneLineSummary = "Accumulate flag statistics given a BAM file",
    programGroup = DiagnosticsAndQCProgramGroup.class
)
public final class FlagStat extends ReadWalker {

    private final FlagStatus sum = new FlagStatus();

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        sum.add(read);
    }

    @Override
    public Object onTraversalSuccess() {
        return sum;
    }

    // what comes out of the flagstat
    public static final class FlagStatus implements Serializable {
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

        public FlagStatus merge(final FlagStatus that){
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
            return this;
        }

        public FlagStatus add( final GATKRead read ) {
            this.readCount++;

            if ( read.failsVendorQualityCheck() ) {
                this.QC_failure++;
            }
            if ( read.isDuplicate() ) {
                this.duplicates++;
            }
            if ( ! read.isUnmapped() ) {
                this.mapped++;
            }
            if ( read.isPaired() ) {
                this.paired_in_sequencing++;

                if ( read.isSecondOfPair() ) {
                    this.read2++;
                }
                else if ( read.isFirstOfPair() ) {
                    this.read1++;
                }

                if ( read.isProperlyPaired() ) {
                    this.properly_paired++;
                }

                if ( ! read.isUnmapped() && ! read.mateIsUnmapped() ) {
                    this.with_itself_and_mate_mapped++;

                    if ( ! read.getContig().equals(read.getMateContig()) ) {
                        this.with_mate_mapped_to_a_different_chr++;

                        if ( read.getMappingQuality() >= 5 ) {
                            this.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5++;
                        }
                    }
                }

                if ( ! read.isUnmapped() && read.mateIsUnmapped() ) {
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
            return with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 == that.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5;

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

