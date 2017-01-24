package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.filter.OverclippedReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.function.Predicate;

public class LinkedReadAnalysisFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;

    Predicate<GATKRead> filter;

    public LinkedReadAnalysisFilter(final double minEntropy) {
        final Predicate<GATKRead> readPredicate = ReadFilterLibrary.MAPPED
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK)
                .and(ReadFilterLibrary.NOT_DUPLICATE)
                .and(ReadFilterLibrary.PRIMARY_ALIGNMENT)
                .and(read -> !read.isSupplementaryAlignment())
                .and(read -> read.hasAttribute("BX"))
                .and(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO)
                .and(read -> !filterOut(read, 36, true));
        this.filter = readPredicate;
        if (minEntropy > 0) {
            this.filter = this.filter.and(new ReadEntropyFilter(minEntropy));
        }
    }


    // Stolen from htsjdk.samtools.filter.OverclippedReadFilter so that I can use GATKRead
    public boolean filterOut(final GATKRead record, final int unclippedBasesThreshold, final boolean filterSingleEndClips) {
        int alignedLength = 0;
        int softClipBlocks = 0;
        int minSoftClipBlocks = filterSingleEndClips ? 1 : 2;
        CigarOperator lastOperator = null;

        for ( final CigarElement element : record.getCigar().getCigarElements() ) {
            if ( element.getOperator() == CigarOperator.S ) {
                //Treat consecutive S blocks as a single one
                if(lastOperator != CigarOperator.S){
                    softClipBlocks += 1;
                }

            } else if ( element.getOperator().consumesReadBases() ) {   // M, I, X, and EQ (S was already accounted for above)
                alignedLength += element.getLength();
            }
            lastOperator = element.getOperator();
        }

        return(alignedLength < unclippedBasesThreshold && softClipBlocks >= minSoftClipBlocks);
    }


    @Override
    public boolean test(final GATKRead read) {
        return filter.test(read);
    }
}
