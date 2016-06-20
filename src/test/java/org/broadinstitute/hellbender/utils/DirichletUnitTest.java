package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.util.MathArrays;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 * Created by davidben on 7/25/16.
 */
public class DirichletUnitTest {
    @Test
    public void testSymmetricDirichlet() {
        for (final int numStates : new int[] {1,3,10}) {
            for (final double concentration : new double[] {0.1, 1.0, 10.0}) {
                final Dirichlet dirichlet = Dirichlet.symmetricDirichlet(numStates, concentration);
                Assert.assertEquals(dirichlet.size(), numStates);
                final double[] effectiveWeights = dirichlet.effectiveMultinomialWeights();
                Arrays.stream(dirichlet.effectiveMultinomialWeights())
                        .forEach(x -> Assert.assertEquals(x, effectiveWeights[0], 1e-8));
            }
        }
    }

    @Test
    public void testAdditivePseudocounts() {
        final double[] array1 = {0.1, 2.0, 0.5, 3.0, 1.0};
        final double[] array2 = {1.0, 3.0, 0.5, 0.7, 0.9};

        final Dirichlet dirichlet1 = new Dirichlet(new Dirichlet(array1), array2);
        final Dirichlet dirichlet2 = new Dirichlet(new Dirichlet(array2), array1);
        final Dirichlet dirichlet3 = new Dirichlet(MathArrays.ebeAdd(array1, array2));

        Assert.assertEquals(dirichlet1.effectiveMultinomialWeights(), dirichlet2.effectiveMultinomialWeights());
        Assert.assertEquals(dirichlet1.effectiveMultinomialWeights(), dirichlet3.effectiveMultinomialWeights());
    }

    // if params are large, Dirichlet concentrates about a particular multinomial and thus effective weights converge toward the mode
    @Test
    public void testAsymptoticEffectiveWeights() {
        final List<double[]> parameterArrays = Arrays.asList(new double[] {100, 200, 300}, new double[] {1000, 1, 2, 4});
        for (final double[] params : parameterArrays) {
            final double[] normalizedWeights = MathUtils.normalizeFromRealSpace(new Dirichlet(params).effectiveMultinomialWeights());
            final double sum = MathUtils.sum(params);
            final double[] normalizedParams = Arrays.stream(params).map(x -> x / sum).toArray();
            IntStream.range(0, params.length).forEach(n -> Assert.assertEquals(normalizedParams[n], normalizedWeights[n], 0.01));
        }
    }
}