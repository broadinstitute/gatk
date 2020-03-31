package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.Cigar;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.Set;

/**
 * Filters out reads above a threshold identity (number of matches less deletions), given in bases. Reads mapped to
 * any of the excluded contigs, if provided, will automatically pass the filter regardless of mapping identity.
 */
public final class HostAlignmentReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1L;

    private final int minIdentity;
    private final Set<String> excludedContigs;

    public HostAlignmentReadFilter(final int minIdent, final Set<String> excludedContigs) {
        this.minIdentity = minIdent;
        this.excludedContigs = excludedContigs;
    }

    public HostAlignmentReadFilter(final int minIdent) {
        this.minIdentity = minIdent;
        this.excludedContigs = Collections.emptySet();
    }

    @Override
    public boolean test(final GATKRead read) {
        if ( read.isUnmapped() ) return true;
        if ( excludedContigs.contains(read.getContig()) ) return true;
        final Integer nmVal = read.getAttributeAsInteger("NM");
        if ( nmVal == null ) return true;
        return test(read.getCigar(), nmVal.intValue());
    }

    public boolean test(final Cigar cigar, final int numMismatches) {
        return PSUtils.getMatchesLessDeletions(cigar, numMismatches) < minIdentity;
    }

}
