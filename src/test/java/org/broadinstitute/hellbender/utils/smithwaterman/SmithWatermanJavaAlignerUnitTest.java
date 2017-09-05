package org.broadinstitute.hellbender.utils.smithwaterman;

public final class SmithWatermanJavaAlignerUnitTest extends SmithWatermanAlignerAbstractUnitTest {

    @Override
    protected SmithWatermanJavaAligner getAligner() {
        return SmithWatermanJavaAligner.getInstance();
    }

}
