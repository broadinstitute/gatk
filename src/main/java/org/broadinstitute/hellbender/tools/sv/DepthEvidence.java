package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/** Read counts for an indefinite number of samples on some interval. */
public final class DepthEvidence implements SVFeature {
    private final String contig;
    private final int start;
    private final int end;
    private final int[] counts;

    public static final String BCI_VERSION = "1.0";
    public static final int MISSING_DATA = -1;

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
    public DepthEvidence extractSamples( final Set<String> sampleNames, final Object headerObj ) {
        if ( !(headerObj instanceof SVFeaturesHeader) ) {
            throw new UserException("DepthEvidence feature source without a header.  " +
                                    "We don't know which samples we have.");
        }
        final SVFeaturesHeader header = (SVFeaturesHeader)headerObj;
        final int nCounts = sampleNames.size();
        final int[] newCounts = new int[nCounts];
        int idx = 0;
        for ( final String sampleName : sampleNames ) {
            final Integer sampleIndex = header.getSampleIndex(sampleName);
            newCounts[idx++] = sampleIndex == null ? MISSING_DATA : counts[sampleIndex];
        }
        return new DepthEvidence(contig, start, end, newCounts);
    }

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

    @Override public String toString() {
        final StringBuilder sb = new StringBuilder(contig + "\t" + start + "\t" + end);
        for ( final int count : counts ) {
            sb.append("\t");
            sb.append(count);
        }
        return  sb.toString();
    }
}
