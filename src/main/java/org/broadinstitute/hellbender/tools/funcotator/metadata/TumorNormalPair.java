package org.broadinstitute.hellbender.tools.funcotator.metadata;

import org.apache.commons.lang3.tuple.Pair;

/**
 * Convenience class for tumor normal pair.  Only stores the sample names.
 */
public class TumorNormalPair {

    final private Pair<String, String> pair;

    public TumorNormalPair(final String tumorSampleName, final String normalSampleName) {
        this.pair = Pair.of(tumorSampleName, normalSampleName);
    }

    public String getTumor() {
        return pair.getLeft();
    }

    public String getNormal() {
        return pair.getRight();
    }

    @Override
    public String toString() {
        return "Tumor: " + getTumor() + ", Normal: " + getNormal();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        TumorNormalPair that = (TumorNormalPair) o;

        return pair != null ? pair.equals(that.pair) : that.pair == null;
    }

    @Override
    public int hashCode() {
        return pair != null ? pair.hashCode() : 0;
    }
}
