package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Assert;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

/**
 * Created by David Benjamin on 3/9/17.
 */
public class SomaticLikelihoodsEngineUnitTest extends BaseTest {
    @Test
    public void testAlleleFractionsPosterior() {

        //likelihoods completely favor allele 0 over allele 1 for every read, so
        // we should get no counts for allele 1
        final double[] prior1 = new double[] {1, 1};
        final RealMatrix mat1 = new Array2DRowRealMatrix(2, 4);
        mat1.setRow(0, new double[] {0, 0, 0, 0});
        mat1.setRow(1, new double[] {-10, -10, -10, -10});
        final double[] posterior1 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat1, prior1);
        final double[] expectedCounts1 = new double[] {4, 0};

        final double expectedPosterior1[] = MathArrays.ebeAdd(prior1, expectedCounts1);
        Assert.assertArrayEquals(posterior1, expectedPosterior1, 1.0e-6);

        //prior is extremely strong and outweighs ambiguous likelihoods
        final double[] prior2 = new double[] {1e8, 1};
        final RealMatrix mat2 = new Array2DRowRealMatrix(2, 4);
        mat2.setRow(0, new double[] {0, 0, 0, 0});
        mat2.setRow(1, new double[] {0, 0, 0, 0});
        final double[] posterior2 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat2, prior2);
        final double[] expectedCounts2 = new double[] {4, 0};

        final double expectedPosterior2[] = MathArrays.ebeAdd(prior2, expectedCounts2);
        Assert.assertArrayEquals(posterior2, expectedPosterior2, 1.0e-6);

        //prior is extremely weak and likelihoods speak for themselves
        final double[] prior3 = new double[] {1e-6, 1e-6};
        final RealMatrix mat3 = new Array2DRowRealMatrix(2, 4);
        mat3.setRow(0, new double[] {0, 0, 0, -10});
        mat3.setRow(1, new double[] {-10, -10, -10, 0});
        final double[] posterior3 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat3, prior3);
        final double[] expectedCounts3 = new double[] {3, 1};

        final double expectedPosterior3[] = MathArrays.ebeAdd(prior3, expectedCounts3);
        Assert.assertArrayEquals(posterior3, expectedPosterior3, 1.0e-6);

        // test convergence
        final double[] prior4 = new double[] {0.2, 1.7};
        final RealMatrix mat4 = new Array2DRowRealMatrix(2, 4);
        mat4.setRow(0, new double[] {0.1, 5.2, 0.5, 0.2});
        mat4.setRow(1, new double[] {2.6, 0.6, 0.5, 0.4});
        final double[] posterior4 = SomaticLikelihoodsEngine.alleleFractionsPosterior(mat4, prior4);
        final double[] counts4 = MathArrays.ebeSubtract(posterior4, prior4);
        Assert.assertArrayEquals(counts4, SomaticLikelihoodsEngine.getEffectiveCounts(mat4, posterior4), 1.0e-3);
    }

    @Test
    public void testDirichletNormalization() {
        // a Dirichlet with two parameters is a Beta
        final double a = 11;
        final double b = 36;
        final double[] params1 = new double[] {a, b};
        final double normalization = Math.pow(10, SomaticLikelihoodsEngine.log10DirichletNormalization(params1));
        final BetaDistribution bd = new BetaDistribution(a,b);
        for (final double x : new double[] {0.1, 0.3, 0.6}) {
            Assert.assertEquals(bd.density(x), normalization * Math.pow(x, a - 1) * Math.pow(1-x, b-1), 1e-6);
        }

        // a Dirichlet with parameters equal to 1 is flat, so the normalization is 1/the hypervolume of the simplex
        // which is d! where d is the dimension

        final double[] params2 = new double[] {1, 1, 1, 1};
        final double normalization2 = Math.pow(10, SomaticLikelihoodsEngine.log10DirichletNormalization(params2));
        Assert.assertEquals(normalization2, 6, 1e-6);
    }

    @Test
    public void testEvidence() {
        // one exact limit for the evidence is when the likelihoods of each read are so peaked (i.e. the most likely allele
        // of each read is much likelier than all other alleles) that the sum over latent read-to-allele assignments
        // (that is, over the indicator z in the notes) is dominated by the max-likelihood allele configuration
        // and thus the evidence reduces to exactly integrating out the Dirichlet allele fractions

        final double[] prior = new double[] {1, 2};
        final RealMatrix log10Likelihoods = new Array2DRowRealMatrix(2, 4);
        log10Likelihoods.setRow(0, new double[] {0.1, 4.0, 3.0, -10});
        log10Likelihoods.setRow(1, new double[] {-12, -9, -5.0, 0.5});
        final double calculatedLog10Evidence = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods, prior);
        final double[] maxLikelihoodCounts = new double[] {3, 1};
        final double expectedLog10Evidence = SomaticLikelihoodsEngine.log10DirichletNormalization(prior)
                - SomaticLikelihoodsEngine.log10DirichletNormalization(MathArrays.ebeAdd(prior, maxLikelihoodCounts))
                + new IndexRange(0,log10Likelihoods.getColumnDimension()).sum(read -> log10Likelihoods.getColumnVector(read).getMaxValue());
        Assert.assertEquals(calculatedLog10Evidence, expectedLog10Evidence, 1e-5);

        // when there's just one read we can calculate the likelihood exactly

        final double[] prior2 = new double[] {1, 2};
        final RealMatrix log10Likelihoods2 = new Array2DRowRealMatrix(2, 1);
        log10Likelihoods2.setRow(0, new double[] {0.1});
        log10Likelihoods2.setRow(1, new double[] {0.5});
        final double calculatedLog10Evidence2 = SomaticLikelihoodsEngine.log10Evidence(log10Likelihoods2, prior2);
        final double[] delta0 = new double[] {1, 0};
        final double[] delta1 = new double[] {0, 1};
        final double expectedLog10Evidence2 = MathUtils.log10SumLog10(log10Likelihoods2.getEntry(0,0) +
                SomaticLikelihoodsEngine.log10DirichletNormalization(prior2)
                - SomaticLikelihoodsEngine.log10DirichletNormalization(MathArrays.ebeAdd(prior2, delta0)),
                + log10Likelihoods2.getEntry(1,0) +
                SomaticLikelihoodsEngine.log10DirichletNormalization(prior2)
                        - SomaticLikelihoodsEngine.log10DirichletNormalization(MathArrays.ebeAdd(prior2, delta1)));
        Assert.assertEquals(calculatedLog10Evidence2, expectedLog10Evidence2, 0.05);


    }

}