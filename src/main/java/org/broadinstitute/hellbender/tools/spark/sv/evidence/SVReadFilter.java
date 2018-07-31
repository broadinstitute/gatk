package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Iterator;
import java.util.function.BiPredicate;

public class SVReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    private final int minEvidenceMapQ;
    private final int minEvidenceMatchLength;
    private final int allowedShortFragmentOverhang;
    private final int maxDUSTScore;
    private final int kSize;

    public SVReadFilter( final FindBreakpointEvidenceSparkArgumentCollection params ) {
        minEvidenceMapQ = params.minEvidenceMapQ;
        minEvidenceMatchLength = params.minEvidenceMatchLength;
        allowedShortFragmentOverhang = params.allowedShortFragmentOverhang;
        maxDUSTScore = params.maxDUSTScore;
        kSize = params.kSize;
    }

    public boolean notJunk( final GATKRead read ) {
        return !read.isDuplicate() && !read.failsVendorQualityCheck();
    }

    public boolean isPrimaryLine( final GATKRead read ) {
        return !read.isSecondaryAlignment() && !read.isSupplementaryAlignment();
    }

    public boolean isMapped( final GATKRead read ) {
        return notJunk(read) && !read.isUnmapped();
    }

    public boolean isMappedPrimary( final GATKRead read ) {
        return isMapped(read) && isPrimaryLine(read);
    }

    public boolean isEvidence(final GATKRead read) {
        return isMapped(read) &&
                read.getMappingQuality() >= minEvidenceMapQ &&
                CigarUtils.countAlignedBases(read.getCigar()) >= minEvidenceMatchLength && ! read.isSecondaryAlignment();
    }

    public boolean isMappedToPrimaryContig(final GATKRead read, final ReadMetadata readMetadata) {
        if (! isMapped(read)) return false;
        final int contigID = readMetadata.getContigID(read.getContig());
        return ( ! readMetadata.getCrossContigIgnoreSet().contains(contigID));
    }

    public boolean isTemplateLenTestable( final GATKRead read ) {
        return isEvidence(read) && isPrimaryLine(read) &&
                !read.mateIsUnmapped() &&
                !read.isReverseStrand() &&
                read.mateIsReverseStrand() &&
                read.getContig().equals(read.getMateContig()) &&
                read.getStart() - allowedShortFragmentOverhang <= read.getMateStart();
    }

    public boolean containedInRegionToIgnore(final SVInterval interval, final SVIntervalTree<SVInterval> regionsToIgnore) {
        final Iterator<SVIntervalTree.Entry<SVInterval>> overlappers = regionsToIgnore.overlappers(interval);
        while (overlappers.hasNext()) {
            SVIntervalTree.Entry<SVInterval> depthFilteredInterval = overlappers.next();
            if (depthFilteredInterval.getInterval().overlapLen(interval) == interval.getLength()) {
                return true;
            }
        }
        return false;
    }

    public Iterator<GATKRead> applyFilter( final Iterator<GATKRead> readItr, final BiPredicate<SVReadFilter, GATKRead> predicate ) {
        return new SVUtils.IteratorFilter<>(readItr, read -> predicate.test(this, read));
    }

    public int getMinEvidenceMapQ() {
        return minEvidenceMapQ;
    }

    public int getMinEvidenceMatchLength() {
        return minEvidenceMatchLength;
    }

    public int getAllowedShortFragmentOverhang() {
        return allowedShortFragmentOverhang;
    }
}
