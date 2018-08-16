package org.broadinstitute.hellbender.testutils;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.nio.charset.Charset;

public class FakeReferenceSource {
    /**
     * bases returns some made up bases for an interval. The pattern is
     * tagctagctagctagctagc...
     * To make the fake bases consistant, we offset the pattern based on the interval. For example, assuming
     * the contig is 0-based, the interval [2,8) would be
     * gctagc
     * @param interval the interval on the fake contig
     * @return the fake bases.
     */
    public static ReferenceBases bases(SimpleInterval interval) {
        int start = interval.getStart();
        int end = interval.getEnd();

        int chunkStart = 4*(Math.floorDiv(start, 4));
        int chunkEnd = 4*(Math.floorDiv(end + 1, 4) + 1);
        String full = StringUtils.repeat("tagc", chunkEnd - chunkStart);
        String substring = full.substring(start - chunkStart, end - chunkStart + 1);

        return new ReferenceBases(substring.getBytes(Charset.forName("UTF-8")), interval);
    }

}
