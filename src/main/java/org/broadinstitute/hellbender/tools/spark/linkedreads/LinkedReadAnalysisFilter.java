package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class LinkedReadAnalysisFilter implements ReadFilter {
    private static final long serialVersionUID = 1l;

    ReadFilter filter;

    public LinkedReadAnalysisFilter() {
        this.filter = ReadFilterLibrary.MAPPED
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK)
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK)
                .and(GATKRead::isFirstOfPair)
                .and(ReadFilterLibrary.NOT_DUPLICATE)
                .and(ReadFilterLibrary.PRIMARY_ALIGNMENT)
                .and(read -> !read.isSupplementaryAlignment())
                .and(read -> read.hasAttribute("BX"));
    }

    @Override
    public boolean test(final GATKRead read) {
        return filter.test(read);
    }
}
