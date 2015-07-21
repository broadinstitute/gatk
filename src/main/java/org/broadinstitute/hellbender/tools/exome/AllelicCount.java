package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Interval;
import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * Reference and alternate allele counts at a SNP site specified by an interval.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCount {
    private final Interval interval;
    private final int refReadCount, altReadCount;

    public AllelicCount(final Interval interval, final int refReadCount, final int altReadCount) {
        this.interval = interval;
        this.refReadCount = refReadCount;
        this.altReadCount = altReadCount;
    }

    public Interval getInterval() { return interval.clone();    }

    public int getRefReadCount() {  return refReadCount;        }

    public int getAltReadCount() {  return altReadCount;        }

    public TargetCoverage asTargetCoverage(final String name, final double allelicFractionSkew) {
        final double coverage = Math.abs((double) altReadCount / (refReadCount + altReadCount) - allelicFractionSkew/2);
        final TargetCoverage countAsTargetCoverage = new TargetCoverage(name, new SimpleInterval(interval),
                coverage);
        return countAsTargetCoverage;
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
