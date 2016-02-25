package org.broadinstitute.hellbender.tools.spark.sv;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Customized serializable reads filter, based on cmd line arguments provided, for structural variation purposes
 */
public final class SVCustomReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;

    private static final PrimaryFilter pfilter = new PrimaryFilter();
    private final CustomFilter sfilter;

    public SVCustomReadFilter(final boolean[] sixFilterOptions,
                              final whichEndToUse whichEnd,
                              final int mappingQualityThreshold
                              //TODO: final int mateMappingQualityThreshold
                             ){
        sfilter = new CustomFilter(sixFilterOptions, whichEnd, mappingQualityThreshold);//TODO:, mateMappingQualityThreshold;
    }
    @Override
    public boolean test(final GATKRead read){
        return pfilter.and(sfilter).test(read);
    }

    /**
     * Customary reads filter where filters are turned on/off by user options, which are passed in via constructor.
     */
    private static final class CustomFilter implements ReadFilter{
        private static final long serialVersionUID = 1L;

        private final int endVal;

        private final ReadFilter collapsedFilter;

        //TODO: private final int mateMQPassingThreshold = 0;

        public CustomFilter(final boolean sixOptions[], final whichEndToUse whichEnd, final int MQThreshold//, mateMQThreshold
                            ){

            endVal = whichEnd.getValue();

            ReadFilter tempFilter = read -> 0!=read.getFragmentLength();
                       tempFilter = tempFilter.and(read -> read.isPaired());
                       tempFilter = tempFilter.and(read -> endVal == (read.isFirstOfPair() ? 1 : 2));

            if(sixOptions[0]) { tempFilter = tempFilter.and(read -> read.isProperlyPaired());}
            if(sixOptions[1]) { tempFilter = tempFilter.and(read -> !read.isUnmapped());}
            if(sixOptions[2]) { tempFilter = tempFilter.and(read -> !read.mateIsUnmapped());}
            if(sixOptions[3]) { tempFilter = tempFilter.and(read -> !read.isDuplicate());}
            if(sixOptions[4]) { tempFilter = tempFilter.and(read -> !read.isSecondaryAlignment());}
            if(sixOptions[5]) { tempFilter = tempFilter.and(read -> !read.isSupplementaryAlignment());}
            if(0!=MQThreshold){ tempFilter = tempFilter.and(read -> read.getMappingQuality() >= MQThreshold);}

            collapsedFilter = tempFilter;

            // TODO: efficiently checking mate mapping quality (may require index information)
            // TODO: pick only isReverseStrand/mateIsReverseStrand (strand bias)
            // TODO: case 0 of "whichEnd", most complicated since behavior depends on if mate is available in interesting region
            // TODO: filter or not based on length==0
            // TODO: filter based on MATE_ON_SAME_CONTIG MATE_DIFFERENT_STRAND GOOD_CIGAR NON_ZERO_REFERENCE_LENGTH_ALIGNMENT
        }

        @Override
        public boolean test(final GATKRead read){
            return collapsedFilter.test(read);
        }
    }
}

/**
 * Similar to WellformedReaFilter except doesn't check for header but checks if mapping quality is available.
 */
final class PrimaryFilter implements ReadFilter{
    private static final long serialVersionUID = 1L;

    private final ReadFilter primary = ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK
                                       .and(ReadFilterLibrary.VALID_ALIGNMENT_START)
                                       .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                                       .and(ReadFilterLibrary.HAS_READ_GROUP)
                                       .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                                       .and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                                       .and(ReadFilterLibrary.SEQ_IS_STORED)
                                       .and(ReadFilterLibrary.CIGAR_IS_SUPPORTED)
                                       .and(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);

    @Override
    public boolean test(final GATKRead read){
        return primary.test(read);
    }
}
