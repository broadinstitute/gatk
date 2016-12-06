package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class LinkedReadAnalysisFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;

    public static final int MAX_FRAGMENT_LENGTH = 10000;
    ReadFilter filter;

    public LinkedReadAnalysisFilter() {
        this.filter = ReadFilterLibrary.MAPPED
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK)
                .and(ReadFilterLibrary.NOT_DUPLICATE)
                .and(ReadFilterLibrary.PRIMARY_ALIGNMENT);
    }

    @Override
    public boolean test(final GATKRead read) {
        return filter.test(read) &&
                read.isFirstOfPair() &&
                read.hasAttribute("BX") &&
                !isChimeric(read);
    }

    public static boolean isChimeric(final GATKRead read) {
        return ((read.getMateContig() != null && !read.getContig().equals(read.getMateContig())) ||
                Math.abs(read.getFragmentLength()) >= MAX_FRAGMENT_LENGTH ||
                read.isReverseStrand() == read.mateIsReverseStrand() ||
                (read.getStart() < read.getMateStart() && read.isReverseStrand()) ||
                (read.getStart() > read.getMateStart() && !read.isReverseStrand()));
    }


}
