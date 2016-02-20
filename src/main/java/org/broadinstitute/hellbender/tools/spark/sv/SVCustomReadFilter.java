package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Customized serializable reads filter, based on cmd line arguments provided, for structural variation purposes
 */
public final class SVCustomReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;

    // If set to true, filter out pairs of reads that are not properly oriented. Default value: true.
    private boolean properlyPairedOnly;

    // If set to true, include mapped reads only.
    private boolean mappedReadsOnly;

    // If set to true, include mate-mapped reads only.
    private boolean mateMappedReadsOnly;

    // If set to true, include non-duplicated reads only.
    private boolean nonDupReadsOnly;

    // If set to true, filter out secondary alignments reads.
    private boolean nonSecondaryAlignmentsOnly;

    // If set to true, filter out supplementary alignments reads.
    private boolean nonSupplementaryAlignmentsOnly;

    // If set to true, collect from reads that are first end of a pair.
    private int whichEnd;

    // If set non-zero value, only include reads passing certain mapping quality threshold.
    private double MQPassingThreshold = 0.0;

    // If set to non-zero value, only include reads with mate passing certain mapping quality threshold.
    private double mateMQPassingThreshold = 0.0;

    private final primaryFilter pfilter = new primaryFilter();
    private final customFilter sfilter = new customFilter();


    public SVCustomReadFilter(final boolean[] sixFilterOptions,
                              final int whichEnd,
                              final double mappingQualityThreshold,
                              final double mateMappingQualityThreshold){
        properlyPairedOnly = sixFilterOptions[0];
        mappedReadsOnly = sixFilterOptions[1];
        mateMappedReadsOnly = sixFilterOptions[2];
        nonDupReadsOnly = sixFilterOptions[3];
        nonSecondaryAlignmentsOnly = sixFilterOptions[4];
        nonSupplementaryAlignmentsOnly = sixFilterOptions[5];

        this.whichEnd = whichEnd;

        MQPassingThreshold = mappingQualityThreshold;
        mateMQPassingThreshold = mateMappingQualityThreshold;
    }
    @Override
    public boolean test(final GATKRead read){
        return pfilter.and(sfilter).test(read);
    }

    // similar to WellformedReadFilter except doesn't check for header but checks if mapping quality is available
    private final class primaryFilter implements ReadFilter{
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test(final GATKRead read){
            final ReadFilter primary = ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK
                                        .and(ReadFilterLibrary.VALID_ALIGNMENT_START)
                                        .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                                        .and(ReadFilterLibrary.HAS_READ_GROUP)
                                        .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                                        .and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                                        .and(ReadFilterLibrary.SEQ_IS_STORED)
                                        .and(ReadFilterLibrary.CIGAR_IS_SUPPORTED)
                                        .and(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
            return primary.test(read);
        }
    }

    private final class customFilter implements ReadFilter{
        private static final long serialVersionUID = 1L;

        @Override
        public boolean test(final GATKRead read){

            boolean result = true;

            if(properlyPairedOnly)
                result &= read.isProperlyPaired();
            if(mappedReadsOnly)
                result &= !read.isUnmapped();
            if(mateMappedReadsOnly)
                result &= !read.mateIsUnmapped();
            if(nonDupReadsOnly)
                result &= !read.isDuplicate();
            if(nonSecondaryAlignmentsOnly)
                result &= !read.isSecondaryAlignment();
            if(nonSupplementaryAlignmentsOnly)
                result &= !read.isSupplementaryAlignment();

            if(1==whichEnd)
                result &= read.isFirstOfPair();
            else if(2==whichEnd)
                result &= read.isSecondOfPair();

            if(MQPassingThreshold > 0.0)
                result &= (read.getMappingQuality() >= MQPassingThreshold);
            // TODO: efficiently checking mate mapping quality
            // TODO: pick only isReverseStrand/mateIsReverseStrand (strand bias)
            // TODO: case 0 of "whichEnd", most complicated since behavior depends on if mate is available in interesting region
            // TODO: filter or not based on length==0
            // TODO: filter based on MATE_ON_SAME_CONTIG MATE_DIFFERENT_STRAND GOOD_CIGAR NON_ZERO_REFERENCE_LENGTH_ALIGNMENT
            return result;
        }
    }
}