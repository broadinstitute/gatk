package org.broadinstitute.hellbender.tools.picard.illumina.parser;

import static java.lang.Math.pow;

/**
 * While structurally identical to CompositeIndex, this class is maintained as it makes code more readable when the two are used together (see QSeqParser)
 *
 * @author jburke@broadinstitute.org
 */
public final class Range {
    public final int start;
    public final int end;
    public final int length;

    public Range(final int start, final int end) {
        if (end < start) {
            throw new IlluminaParserException("Nonsensical Range(" + start + ", " + end + ")");
        }

        this.start = start;
        this.end = end;
        this.length = end - start + 1;
    }

    @Override
    public boolean equals(final Object object) {
        if (object == null || !(object instanceof Range)) {
            return false;
        }

        final Range that = (Range) object;
        return that.start == this.start && that.end == this.end;
    }

    @Override
    public int hashCode() {
        return (int) pow(start, end);
    }

    @Override
    public String toString() {
        return "Range(" + start + ", " + end + ", " + length + ")";
    }
}
