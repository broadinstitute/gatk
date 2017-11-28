package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

/**
 * Reference and alternate allele counts at a site specified by an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class AllelicCount implements Locatable {

    /* these are mandatory */
    private final SimpleInterval interval;
    private final int refReadCount;
    private final int altReadCount;
    /* these are set to N if not specified */
    private final Nucleotide refNucleotide;
    private final Nucleotide altNucleotide;

    /**
     * Construct the allelic count object.
     * @param interval the genomic interval (it is assumed to have the same begin and end points)
     * @param refReadCount ref nucleotide count
     * @param altReadCount alt nucleotide count
     * @param refNucleotide ref nucleotide
     * @param altNucleotide alt nucleotide
     */
    public AllelicCount(final SimpleInterval interval,
                        final int refReadCount,
                        final int altReadCount,
                        final Nucleotide refNucleotide,
                        final Nucleotide altNucleotide) {
        this.interval = Utils.nonNull(interval);
        this.refReadCount = ParamUtils.isPositiveOrZero(refReadCount, "Can't construct AllelicCount with negative read counts.");
        this.altReadCount = ParamUtils.isPositiveOrZero(altReadCount, "Can't construct AllelicCount with negative read counts.");
        this.refNucleotide = Utils.nonNull(refNucleotide);
        this.altNucleotide = Utils.nonNull(altNucleotide);
    }

    public AllelicCount(final SimpleInterval interval, final int refReadCount, final int altReadCount) {
        this(interval, refReadCount, altReadCount, Nucleotide.N, Nucleotide.N);
    }

    @Override
    public String getContig() { return interval.getContig(); }

    @Override
    public int getStart() { return interval.getStart(); }

    @Override
    public int getEnd() { return interval.getEnd(); }

    public SimpleInterval getInterval() { return interval; }

    public int getRefReadCount() { return refReadCount; }

    public int getAltReadCount() { return altReadCount; }

    public Nucleotide getRefNucleotide() {
        return refNucleotide;
    }

    public Nucleotide getAltNucleotide() {
        return altNucleotide;
    }

    public int getTotalReadCount() {
        return altReadCount + refReadCount;
    }

    public double getAlternateAlleleFraction() {
        return getTotalReadCount() == 0 ? 0. : (double) altReadCount / getTotalReadCount();
    }

    /**
     * If all fields are specified for both {@link AllelicCount} objects, they are all used to check for equality.
     * Otherwise, if either or both counts have both nucleotides unspecified, then
     * only the intervals and the ref and alt counts are used.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCount)) {
            return false;
        }

        final AllelicCount count = (AllelicCount) o;
        return interval.equals(count.interval)
                && refReadCount == count.refReadCount && altReadCount == count.altReadCount
                && ((refNucleotide == Nucleotide.N && altNucleotide == Nucleotide.N) ||
                (count.refNucleotide == Nucleotide.N && count.altNucleotide == Nucleotide.N) ||
                (refNucleotide == count.refNucleotide && altNucleotide == count.altNucleotide));
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + refReadCount;
        result = 31 * result + altReadCount;
        return result;
    }

    @Override
    public String toString() {
        return "AllelicCount{" +
                "interval=" + interval +
                ", refReadCount=" + refReadCount +
                ", altReadCount=" + altReadCount +
                ", refNucleotide=" + refNucleotide +
                ", altNucleotide=" + altNucleotide +
                '}';
    }
}
