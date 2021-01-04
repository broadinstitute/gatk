package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.Objects;
import org.broadinstitute.hellbender.utils.codecs.DepthEvidenceCodec;

import java.util.ArrayList;
import java.util.List;

public final class DepthEvidence implements Feature {

    final String contig;
    final int start;
    final int end;
    final int[] counts;

    public DepthEvidence(final String contig, int start, final int end, final int[] counts) {
        Utils.nonNull(contig);
        Utils.nonNull(counts);
        this.contig = contig;
        this.start = start;
        this.end = end;
        this.counts = counts;
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public int[] getCounts() { return counts; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof DepthEvidence)) return false;
        DepthEvidence that = (DepthEvidence) o;
        return start == that.start &&
                end == that.end &&
                contig.equals(that.contig) &&
                Arrays.equals(counts, that.counts);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(contig, start, end);
        result = 31 * result + Arrays.hashCode(counts);
        return result;
    }

    public String toString() {
        final int numCounts = counts.length;
        final List<String> data = new ArrayList<>(3 + numCounts);
        data.add(contig);
        data.add(Integer.toString(start - 1));
        data.add(Integer.toString(end));
        for (int i = 0; i < numCounts; i++) {
            data.add(Integer.toString(counts[i]));
        }
        return String.join(DepthEvidenceCodec.COL_DELIMITER, data);
    }
}