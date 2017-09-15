package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.Cigar;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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
        return PSUtils.getMatchesLessDeletions(cigar, numMismatches) < minIdentity;
    }

}
