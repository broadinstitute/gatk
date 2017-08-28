package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Filters out reads above a threshold identity (number of matches less deletions), given in bases.
 */
public final class HostAlignmentReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    private final int minIdentity;

    public HostAlignmentReadFilter(final int minIdent) {
        this.minIdentity = minIdent;
    }

    @Override
    public boolean test(final GATKRead read) {
        if ( read.isUnmapped() ) return true;
        final Integer nmVal = read.getAttributeAsInteger("NM");
        if ( nmVal == null ) return true;
        return test(read.getCigar(), nmVal.intValue());
    }

    public boolean test(final Cigar cigar, final int numMismatches) {
        return getMatchesLessDeletions(cigar, numMismatches) < minIdentity;
    }

    static int getMatchesLessDeletions(final Cigar cigar, final int numMismatches) {
        int numMatches = -numMismatches;
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        for (final CigarElement e : cigarElements) {
            if (e.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH)) {
                numMatches += e.getLength();
            } else if (e.getOperator().equals(CigarOperator.DELETION)) {
                numMatches -= e.getLength();
            }
        }
        return numMatches;
    }

}
