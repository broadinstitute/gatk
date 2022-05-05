package org.broadinstitute.hellbender.tools.sv;


import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;
import java.util.Set;

/** Documents evidence of a too-close or too-far-apart read pair. */
public final class DiscordantPairEvidence implements SVFeature {

    final String sample;
    final String startContig;
    final String endContig;
    final int start;
    final int end;
    final boolean startStrand;
    final boolean endStrand;

    public final static String BCI_VERSION = "1.0";

    public DiscordantPairEvidence(final String sample,
                                  final String startContig, final int start, final boolean startStrand,
                                  final String endContig, final int end, final boolean endStrand) {
        Utils.nonNull(sample);
        Utils.nonNull(startContig);
        Utils.nonNull(endContig);
        this.sample = sample;
        this.startContig = startContig;
        this.start = start;
        this.startStrand = startStrand;
        this.endContig = endContig;
        this.end = end;
        this.endStrand = endStrand;
    }

    public String getSample() {
        return sample;
    }

    // For purposes of indexing, we will return the start position
    @Override
    public String getContig() {
        return startContig;
    }

    @Override
    public int getStart() {
        return start;
    }

    public boolean getStartStrand() {
        return startStrand;
    }

    public String getEndContig() {
        return endContig;
    }

    /**
     * Used for sorting/checking evidence file order, where end coordinates should not be compared.
     * @return start position
     */
    @Override
    public int getEnd() {
        return start;
    }

    public int getEndPosition() {
        return end;
    }

    public boolean getEndStrand() {
        return endStrand;
    }

    @Override
    public DiscordantPairEvidence extractSamples( final Set<String> sampleNames,
                                                  final Object header ) {
        return sampleNames.contains(sample) ? this : null;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof DiscordantPairEvidence)) return false;
        DiscordantPairEvidence that = (DiscordantPairEvidence) o;
        return start == that.start &&
                end == that.end &&
                startStrand == that.startStrand &&
                endStrand == that.endStrand &&
                sample.equals(that.sample) &&
                startContig.equals(that.startContig) &&
                endContig.equals(that.endContig);
    }

    @Override
    public int hashCode() {
        return Objects.hash(sample, startContig, endContig, start, end, startStrand, endStrand);
    }

    @Override public String toString() {
        return startContig + "\t" + start + "\t" + end + "\t" + sample + "\t" + endContig + "\t" +
                startStrand + "\t" + endStrand;
    }
}
