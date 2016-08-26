package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.MathArrays;

import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;

/**
 * The Dirichlet distribution is a distribution on multinomial distributions: if pi is a vector of positive multinomial weights
 * such that sum_i pi[i] = 1, the Dirichlet pdf is P(pi) = [prod_i Gamma(alpha[i]) / Gamma(sum_i alpha[i])] * prod_i pi[i]^(alpha[i] - 1)
 *
 * The vector alpha comprises the sufficient statistics for the Dirichlet distribution.
 *
 * Since the Dirichlet is the conjugate prior to the multinomial, if one has a Dirichlet prior with concentration alpha
 * and observes each category i n_i times (assuming categories are drawn from a multinomial distribution pi)
 * the posterior is alpha_i -> alpha_i + n_i
 *
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class Dirichlet {
    final double[] alpha;

    public Dirichlet(final double... alpha) {
        Utils.nonNull(alpha);
        Utils.validateArg(alpha.length >= 1, "Dirichlet parameters must have at least one element");
        Utils.validateArg(MathUtils.allMatch(alpha, x -> x >= 0), "Dirichlet parameters may not be negative");
        Utils.validateArg(MathUtils.allMatch(alpha, Double::isFinite), "Dirichlet parameters must be finite");
        this.alpha = alpha.clone();
    }

    /**
     * Create a symmetric distribution Dir(a/K, a/K, a/K . . .) where K is the number of states and
     * a is the concentration.
     */
    public static Dirichlet symmetricDirichlet(final int numStates, final double concentration) {
        Utils.validateArg(numStates > 0, "Must have at leat one state");
        Utils.validateArg(concentration > 0, "concentration must be positive");
        return new Dirichlet(Collections.nCopies(numStates, concentration/numStates).stream().mapToDouble(x->x).toArray());
    }

    // in variational Bayes one often needs the effective point estimate of a multinomial distribution with a
    // Dirichlet prior.  This value is not the mode or mean of the Dirichlet but rather the exp of the expected log weights.
    // note that these effective weights do not add up to 1.  This is fine because in any probabilistic model scaling all weights
    // amounts to an arbitrary normalization constant, but it's important to keep in mind because some classes may expect
    // normalized weights.  In that case the calling code must normalize the weights.
    public double[] effectiveMultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return MathUtils.applyToArray(alpha, a -> Math.exp(Gamma.digamma(a) - digammaOfSum));
    }

    public double[] effectiveLog10MultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return MathUtils.applyToArray(alpha, a -> (Gamma.digamma(a) - digammaOfSum) * MathUtils.LOG10_OF_E);
    }

    public double[] meanWeights() {
        final double sum = MathUtils.sum(alpha);
        return MathUtils.applyToArray(alpha, x -> x / sum);
    }

    public double[] log10MeanWeights() {
        final double sum = MathUtils.sum(alpha);
        return MathUtils.applyToArray(alpha, x -> Math.log10(x / sum));
    }

    public int size() { return alpha.length; }
}
