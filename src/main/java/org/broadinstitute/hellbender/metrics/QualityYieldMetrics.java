package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;

/** A set of metrics used to describe the general quality of a BAM file */
public final class QualityYieldMetrics extends MetricBase implements Serializable {
    private static final long serialVersionUID = 1;

    private boolean useOriginalQualities;
    private boolean finalized;

    //Note: this class must have a 0-arg constructor so we need to use a setter for this field.
    //Also, we want to be able to easily construct usable instances so the setter returns 'this' so we can say
    //new QualityYieldMetrics().setUseOriginalQualities(true) in the aggregate call
    public QualityYieldMetrics setUseOriginalQualities(final boolean useOriginalQualities){
        this.useOriginalQualities = useOriginalQualities;
        return this;
    }

    //Note: The names of those fields are non-conventional (upper case)
    //but Picard metrics are expected to have upper-case names
    //Note: those fields must be public because code in superclass finds them only if they are public.

    /** The total number of reads in the input file */
    public int TOTAL_READS;

    /** The number of reads that are PF - pass filter */
    public int PF_READS;

    /** The average read length of all the reads (will be fixed for a lane) */
    public int READ_LENGTH;

    /** The total number of bases in all reads */
    public long TOTAL_BASES;

    /** The total number of bases in all PF reads */
    public long PF_BASES;

    /** The number of bases in all reads that achieve quality score 20 or higher */
    public long Q20_BASES;

    /** The number of bases in PF reads that achieve quality score 20 or higher */
    public long PF_Q20_BASES;

    /** The number of bases in all reads that achieve quality score 30 or higher */
    public long Q30_BASES;

    /** The number of bases in PF reads that achieve quality score 30 or higher */
    public long PF_Q30_BASES;

    /** The sum of quality scores of all bases divided by 20 */
    public long Q20_EQUIVALENT_YIELD;

    /** The sum of quality scores of all bases divided by 20 */
    public long PF_Q20_EQUIVALENT_YIELD;

    //adds stats for the read
    public QualityYieldMetrics addRead(final GATKRead read){
        TOTAL_READS++;
        final int length = read.getLength();

        final boolean isPfRead = !read.failsVendorQualityCheck();
        if (isPfRead) {
            PF_READS++;
            PF_BASES += length;
        }

        TOTAL_BASES += length;

        final byte[] quals;
        if (useOriginalQualities) {
            byte[] tmp = ReadUtils.getOriginalBaseQualities(read);
            if (tmp == null) {
                tmp = read.getBaseQualities();
            }
            quals = tmp;
        } else {
            quals = read.getBaseQualities();
        }

        // add up quals, and quals >= 20
        for (final byte qual : quals) {
            Q20_EQUIVALENT_YIELD += qual;
            if (qual >= 20) {
                Q20_BASES++;
            }
            if (qual >= 30) {
                Q30_BASES++;
            }

            if (isPfRead) {
                PF_Q20_EQUIVALENT_YIELD += qual;
                if (qual >= 20) {
                    PF_Q20_BASES++;
                }
                if (qual >= 30) {
                    PF_Q30_BASES++;
                }
            }
        }
        return this;
    }

    //combines two objects
    public QualityYieldMetrics combine(final QualityYieldMetrics that) {
        Utils.nonNull(that);
        if (useOriginalQualities != that.useOriginalQualities){
            throw new IllegalArgumentException("must have the same value for useOriginalQualities");
        }
        if (finalized || that.finalized){
            throw new IllegalArgumentException("can;t merge objects when the calculations are finalized");
        }

        TOTAL_READS  += that.TOTAL_READS;
        PF_READS     += that.PF_READS;
        READ_LENGTH  += that.READ_LENGTH;
        TOTAL_BASES  += that.TOTAL_BASES;
        PF_BASES     += that.PF_BASES;
        Q20_BASES    += that.Q20_BASES;
        PF_Q20_BASES += that.PF_Q20_BASES;
        Q30_BASES    += that.Q30_BASES;
        PF_Q30_BASES += that.PF_Q30_BASES;
        Q20_EQUIVALENT_YIELD += that.Q20_EQUIVALENT_YIELD;
        PF_Q20_EQUIVALENT_YIELD += that.PF_Q20_EQUIVALENT_YIELD;

        return this;
    }

    //completes the calculations
    public QualityYieldMetrics finish() {
        READ_LENGTH = TOTAL_READS == 0 ? 0 : (int) (TOTAL_BASES / TOTAL_READS);
        Q20_EQUIVALENT_YIELD = Q20_EQUIVALENT_YIELD / 20;
        PF_Q20_EQUIVALENT_YIELD = PF_Q20_EQUIVALENT_YIELD / 20;

        finalized = true;
        return this;
    }

    /**
     * Return a disambiguating name suffix to be used by the multiple collectors to
     * decorate output names, which are provided by the user in the form of a base name
     * that needs to be disambiguated for each individual collector. This prevents
     * a collector from clobbering the previous collectors output file(s). The value
     * is not used here.
     */
    public static String getUniqueNameSuffix() { return "qualityYieldMetrics";}

}
