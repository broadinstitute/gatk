package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Reference and alternate allele counts at a SNP site specified by an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCount implements Locatable {
    private final Interval interval;
    private final int refReadCount, altReadCount;

    public AllelicCount(final Interval interval, final int refReadCount, final int altReadCount) {
        Utils.nonNull(interval, "Can't construct AllelicCount with null interval.");
        if (refReadCount < 0 || altReadCount < 0) {
            throw new IllegalArgumentException("Can't construct AllelicCount with negative read counts.");
        }
        this.interval = interval;
        this.refReadCount = refReadCount;
        this.altReadCount = altReadCount;
    }

    @Override
    public String getContig() {return interval.getContig(); }

    @Override
    public int getStart() {return interval.getStart(); }

    @Override
    public int getEnd() {return interval.getEnd(); }

    public Interval getInterval() { return interval;    }

    public int getRefReadCount() {  return refReadCount;        }

    public int getAltReadCount() {  return altReadCount;        }

    public TargetCoverage toTargetCoverage(final String name, final double allelicFractionSkew) {
        final double coverage = Math.abs((double) altReadCount / (refReadCount + altReadCount) - allelicFractionSkew/2);
        final TargetCoverage target = new TargetCoverage(name, new SimpleInterval(interval), coverage);
        return target;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCount)) {
            return false;
        }

        final AllelicCount count = (AllelicCount) o;
        return interval.equals(count.interval) && refReadCount == count.refReadCount
                && altReadCount == count.altReadCount;
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + refReadCount;
        result = 31 * result + altReadCount;
        return result;
    }
}
