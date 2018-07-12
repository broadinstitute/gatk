package org.broadinstitute.hellbender.tools.spark.sv.utils;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

public class FermiLiteCoverageInterpreter {
    private final static int COVERAGE_MASK = 0xff;
    private final static int COVERAGE_OFFSET = 33;

    private final byte[] perBaseCoverage;

    public FermiLiteCoverageInterpreter( final int size ) {
        perBaseCoverage = new byte[size];
        Arrays.fill(perBaseCoverage, (byte)COVERAGE_OFFSET);
    }

    public FermiLiteCoverageInterpreter( final byte[] perBaseCoverage ) {
        Utils.validateArg(perBaseCoverage != null, "perBaseCoverage must not be null");
        this.perBaseCoverage = perBaseCoverage;
    }

    public int size() { return perBaseCoverage.length; }

    public int get( final int idx ) {
        return (perBaseCoverage[idx] & COVERAGE_MASK) - COVERAGE_OFFSET;
    }
    public byte[] getBytes() { return perBaseCoverage; }

    public void set( final int idx, final int value ) {
        perBaseCoverage[idx] = (byte)Math.min(COVERAGE_MASK, value + COVERAGE_OFFSET);
    }

    public void sumCoverage( final FermiLiteCoverageInterpreter source, final int destinationOffset ) {
        int destIdx = destinationOffset;
        final int sourceSize = source.size();
        for ( int idx = 0; idx != sourceSize; ++idx ) {
            set(destIdx, get(destIdx)+source.get(idx));
            destIdx += 1;
        }
    }

    public void sumCoverageReverse( final FermiLiteCoverageInterpreter source, final int destinationOffset ) {
        final int sourceSize = source.size();
        int destIdx = destinationOffset + sourceSize;
        for ( int idx = 0; idx != sourceSize; ++idx ) {
            destIdx -= 1;
            set(destIdx, get(destIdx)+source.get(idx));
        }
    }
}
