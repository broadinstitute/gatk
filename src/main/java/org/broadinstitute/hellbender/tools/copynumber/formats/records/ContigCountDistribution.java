package org.broadinstitute.hellbender.tools.copynumber.formats.records;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * Represents a coverage distribution on a contig as a histogram of the counts up to a maximum count.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class ContigCountDistribution {
    private final String contig;
    private final Map<Integer, Integer> countDistribution;

    public ContigCountDistribution(final String contig,
                                   final Map<Integer, Integer> countDistribution) {
        this.contig = Utils.nonEmpty(contig);
        this.countDistribution = Utils.nonNull(countDistribution);
        Utils.nonEmpty(countDistribution.entrySet());
    }

    public String getContig() {
        return contig;
    }

    public Map<Integer, Integer> getCountDistribution() {
        return countDistribution;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final ContigCountDistribution that = (ContigCountDistribution) o;

        return contig.equals(that.contig) && countDistribution.equals(that.countDistribution);
    }

    @Override
    public int hashCode() {
        int result = contig.hashCode();
        result = 31 * result + countDistribution.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "ContigCountDistribution{" +
                "contig='" + contig + '\'' +
                ", countDistribution=" + countDistribution +
                '}';
    }
}
