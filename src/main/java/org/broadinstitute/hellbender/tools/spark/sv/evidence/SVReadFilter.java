package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Iterator;
import java.util.function.BiPredicate;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;

public class SVReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    private final int minEvidenceMapQ;
    private final int minEvidenceMatchLength;
    private final int allowedShortFragmentOverhang;

    public SVReadFilter( final FindBreakpointEvidenceSparkArgumentCollection params ) {
        minEvidenceMapQ = params.minEvidenceMapQ;
        minEvidenceMatchLength = params.minEvidenceMatchLength;
        allowedShortFragmentOverhang = params.allowedShortFragmentOverhang;
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

    public boolean isEvidence( final GATKRead read ) {
        return isMapped(read) && read.getMappingQuality() >= minEvidenceMapQ &&
                CigarUtils.countAlignedBases(read.getCigar()) >= minEvidenceMatchLength && ! read.isSecondaryAlignment();
    }

    public boolean isTemplateLenTestable( final GATKRead read ) {
        return isEvidence(read) && isPrimaryLine(read) &&
                !read.mateIsUnmapped() &&
                !read.isReverseStrand() &&
                read.mateIsReverseStrand() &&
                read.getContig().equals(read.getMateContig()) &&
                read.getStart() - allowedShortFragmentOverhang <= read.getMateStart();
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
