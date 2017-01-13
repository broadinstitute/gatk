package org.broadinstitute.hellbender.tools.coveragemodel.cachemanager;

import junit.framework.AssertionFailedError;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.testng.annotations.Test;

import java.util.Map;
import java.util.function.Function;

/**
 * Unit tests for {@link ImmutableComputableGraph}
 *
 * TODO github/gatk-protected issue #803
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraphUnitTest extends BaseTest {

    public static Function<Map<String, ? extends Duplicable>, ? extends Duplicable> func_0 = col -> {
        final INDArray x = DuplicableNDArray.of(col.get("X"));
        final double y = DuplicableNumber.of(col.get("Y"));
        return new DuplicableNDArray(x.mul(y));
    };

    public static Function<Map<String, ? extends Duplicable>, ? extends Duplicable> func_1 = col -> {
        final INDArray xProdY = DuplicableNDArray.of(col.get("X_prod_Y"));
        final double y = DuplicableNumber.of(col.get("Y"));
        return new DuplicableNDArray(xProdY.add(y));
    };

    @Test
    public void testCyclicGraph() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testOutdatedCaches() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testUpToDateCaches() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testComputeOnDemandNodes() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testPrimitiveUpdating() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testExternallyComputableUpdating() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testCacheByTag() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testCacheByNode() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testCacheAutoUpdate() {
        throw new AssertionFailedError("Test is not implemented yet");
    }

    @Test
    public void testUnchangedNodesSameReferenceAfterUpdate() {
        throw new AssertionFailedError("Test is not implemented yet");
    }
}
