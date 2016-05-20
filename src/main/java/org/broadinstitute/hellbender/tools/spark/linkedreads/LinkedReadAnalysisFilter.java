package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.filter.OverclippedReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class LinkedReadAnalysisFilter implements ReadFilter {
    private static final long serialVersionUID = 1l;

    ReadFilter filter;

    public LinkedReadAnalysisFilter(final double minEntropy, final SAMFileHeader headerForReads) {
        final OverclippedReadFilter overclippedReadFilter = new OverclippedReadFilter(32, false);
        this.filter = ReadFilterLibrary.MAPPED
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK)
                .and(ReadFilterLibrary.NOT_DUPLICATE)
                .and(ReadFilterLibrary.PRIMARY_ALIGNMENT)
                .and(read -> !read.isSupplementaryAlignment())
                .and(read -> read.hasAttribute("BX"))
                .and(read -> ! overclippedReadFilter.filterOut(read.convertToSAMRecord(headerForReads)));
        if (minEntropy > 0) {
            this.filter = this.filter.and(new ReadEntropyFilter(minEntropy));
        }
    }



    @Override
    public boolean test(final GATKRead read) {
        return filter.test(read);
    }
}
