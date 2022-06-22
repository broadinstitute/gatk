package org.broadinstitute.hellbender.tools.sv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Objects;

public class SVLocus implements SVLocatable {

    private final String contigA;
    private final String contigB;
    private final int positionA;
    private final int positionB;

    public SVLocus(final String contigA, final int positionA, final String contigB, final int positionB) {
        this.contigA = Utils.nonNull(contigA);
        this.contigB = Utils.nonNull(contigB);
        this.positionA = positionA;
        this.positionB = positionB;
    }

    @Override
    public String getContigA() {
        return contigA;
    }

    @Override
    public int getPositionA() {
        return positionA;
    }

    @Override
    public String getContigB() {
        return contigB;
    }

    @Override
    public int getPositionB() { return positionB; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        SVLocus svLocus = (SVLocus) o;
        return positionA == svLocus.positionA && positionB == svLocus.positionB && contigA.equals(svLocus.contigA) && contigB.equals(svLocus.contigB);
    }

    @Override
    public int hashCode() {
        return Objects.hash(contigA, contigB, positionA, positionB);
    }
}
