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
    /* these are extra metadata and can be null */
    private final Nucleotide refNucleotide;
    private final Nucleotide altNucleotide;

    /**
     * Construct the allelic count object.
     * @param interval the genomic interval (it is assumed to have the same begin and end points)
     * @param refReadCount ref nucleotide count in the pileup
     * @param altReadCount alt nucleotide count in the pileup
     * @param refNucleotide ref nucleotide
     * @param altNucleotide alt nucleotide
     */
    public AllelicCount(final SimpleInterval interval, final int refReadCount, final int altReadCount,
                        final Nucleotide refNucleotide, final Nucleotide altNucleotide) {
        this.interval = Utils.nonNull(interval, "Can't construct AllelicCount with null interval.");
        this.refReadCount = ParamUtils.isPositiveOrZero(refReadCount, "Can't construct AllelicCount with negative read counts.");
        this.altReadCount = ParamUtils.isPositiveOrZero(altReadCount, "Can't construct AllelicCount with negative read counts.");
        /* these can be null */
        this.refNucleotide = refNucleotide;
        this.altNucleotide = altNucleotide;
    }

    public AllelicCount(final SimpleInterval interval, final int refReadCount, final int altReadCount) {
        this(interval, refReadCount, altReadCount, null, null);
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
        return Utils.nonNull(refNucleotide, "The ref nucleotide is not set.");
    }

    public Nucleotide getAltNucleotide() {
        return Utils.nonNull(altNucleotide, "The alt nucleotide is not set.");
    }

    /**
     * If all fields are specified for both counts, they are all used to check for equality.
     * Otherwise, only the interval and counts are used.
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
                && (((refNucleotide == null) || (count.refNucleotide == null)) || (refNucleotide == count.refNucleotide))
                && (((altNucleotide == null) || (count.altNucleotide == null)) || (altNucleotide == count.altNucleotide));
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
        return String.format("(%s, r=%s, a=%s)", interval, refReadCount, altReadCount);
    }
}
