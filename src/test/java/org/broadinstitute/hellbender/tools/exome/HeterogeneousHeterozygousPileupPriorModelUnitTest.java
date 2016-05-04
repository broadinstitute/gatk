package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Before;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.stream.IntStream;

/**
 * Unit tests for {@link HeterogeneousHeterozygousPileupPriorModel}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class HeterogeneousHeterozygousPileupPriorModelUnitTest extends BaseTest {

    private static final double MIN_ABNORMAL_FRACTION = 0.1;
    private static final double MAX_ABNORMAL_FRACTION = 0.9;
    private static final int MAX_COPY_NUMBER = 4;
    private static final int QUADRATURE_ORDER = 300;

    private static HeterogeneousHeterozygousPileupPriorModel hetPrior;

    @BeforeClass
    public void initPrior() {
        hetPrior = new HeterogeneousHeterozygousPileupPriorModel(MIN_ABNORMAL_FRACTION, MAX_ABNORMAL_FRACTION,
                        MAX_COPY_NUMBER, QUADRATURE_ORDER);
    }

    @Test
    public void testQuadrature() {
        double sumWeights = hetPrior.gaussIntegrationWeights.stream().mapToDouble(Double::doubleValue).sum();
        double minHetAlleleFraction = (1 - MAX_ABNORMAL_FRACTION) / (MAX_COPY_NUMBER * MAX_ABNORMAL_FRACTION +
                2 * (1 - MAX_ABNORMAL_FRACTION));
        Assert.assertEquals(sumWeights, 1 - 2 * minHetAlleleFraction, 1e-8);
        Assert.assertEquals(hetPrior.gaussIntegrationAbscissas.size(), QUADRATURE_ORDER);
    }

    @Test
    public void testAlelleRatioPrior() {
        /* the prior should integrate to 1 */
        Assert.assertEquals(IntStream.range(0, hetPrior.gaussIntegrationAbscissas.size())
                .mapToDouble(i -> hetPrior.gaussIntegrationWeights.get(i) * hetPrior.alleleFractionPriors.get(i))
                .sum(), 1.0, 1e-3);
        hetPrior.alleleFractionPriors.stream().forEach(x -> Assert.assertTrue(x > 0));
    }

}
