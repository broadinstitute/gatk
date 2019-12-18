package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.Serializable;
import java.util.Collection;

public final class BreakpointPair implements Serializable {
    public static final long serialVersionUID = 1L;
    private final int contigA;
    private final int contigB;
    private final int positionA;
    private final int positionB;
    private final boolean strandA;
    private final boolean strandB;
    private final Collection<String> assembledContigsA;
    private final Collection<String> assembledContigsB;

    public BreakpointPair(final VariantContext mateA, final VariantContext mateB, final SAMSequenceDictionary dictionary) {
        positionA = mateA.getStart();
        positionB = mateB.getStart();
        contigA = dictionary.getSequenceIndex(mateA.getContig());
        contigB = dictionary.getSequenceIndex(mateB.getContig());
        final String altB = mateB.getAlternateAlleles().get(0).toString();
        final String altA = mateA.getAlternateAlleles().get(0).toString();
        if (altB.startsWith("]") || altB.endsWith("]")) {
            strandA = true;
        } else if (altB.startsWith("[") || altB.endsWith("[")) {
            strandA = false;
        } else {
            throw new IllegalArgumentException("Could not get mate strandedness from ALT: " + altB);
        }
        if (altA.startsWith("]") || altA.endsWith("]")) {
            strandB = true;
        } else if (altA.startsWith("[") || altA.endsWith("[")) {
            strandB = false;
        } else {
            throw new IllegalArgumentException("Could not get mate strandedness from ALT: " + altA);
        }
        assembledContigsA = mateA.getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, "");
        assembledContigsB = mateB.getAttributeAsStringList(GATKSVVCFConstants.CONTIG_NAMES, "");
    }

    public Collection<String> getAssembledContigsA() {
        return assembledContigsA;
    }

    public Collection<String> getAssembledContigsB() {
        return assembledContigsB;
    }

    public int getContigA() {
        return contigA;
    }

    public int getContigB() {
        return contigB;
    }

    public int getPositionA() {
        return positionA;
    }

    public int getPositionB() { return positionB; }

    public boolean isStrandA() {
        return strandA;
    }

    public boolean isStrandB() {
        return strandB;
    }

    public PairedStrandedIntervals getPairedStrandedInterval() {
        final StrandedInterval intervalA = new StrandedInterval(new SVInterval(getContigA(), getPositionA(), getPositionA() + 1), isStrandA());
        final StrandedInterval intervalB = new StrandedInterval(new SVInterval(getContigB(), getPositionB(), getPositionB() + 1), isStrandA());
        if (getContigA() < getContigB() || (getContigA() == getContigB() && getPositionA() <= getPositionB())) {
            return new PairedStrandedIntervals(intervalA, intervalB);
        }
        return new PairedStrandedIntervals(intervalB, intervalA);
    }
}
