package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

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
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public class Dirichlet {
    private final double[] alpha;

    /**
     * Create a Dirichlet distribution Dir(alpha[0], alpha[1] . . .)
     */
    public Dirichlet(final double[] alpha) {
        Utils.nonNull(alpha);
        Arrays.stream(alpha).forEach(a -> ParamUtils.isPositive(a, "Dirichlet parameters must be positive."));
        Utils.validateArg(alpha.length > 0, "Dirichlet parameters must have at least one element");
        this.alpha = Arrays.copyOf(alpha, alpha.length);
    }

    /**
     * Create a Dirichlet posterior by summing the Dirichlet prior parameters with observed counts.  This works
     * because the Dirichlet parameters are the sufficient statistics.
     */
    public Dirichlet(final Dirichlet prior, final double[] counts) {
        Utils.nonNull(prior, "prior can't be null");
        Utils.nonNull(counts, "counts can't be null");
        Arrays.stream(counts).forEach(c -> ParamUtils.isPositive(c, "Counts must be positive."));
        Utils.validateArg(counts.length == prior.size(), "Counts and prior must have same length.");
        alpha = IntStream.range(0, prior.size()).mapToDouble(n -> prior.alpha[n] + counts[n]).toArray();
    }

    /**
     * Create a symmetric Dirichlet distribution Dir(a/K, a/K, a/K . . . a/K), where a is the concentration
     * and K is the number of states
     */
    public static Dirichlet symmetricDirichlet(final int numStates, final double concentration) {
        return new Dirichlet(Collections.nCopies(numStates, concentration/numStates).stream().mapToDouble(x -> x).toArray());
    }

    /**
     * in variational Bayes one often needs the effective point estimate of a multinomial distribution with a
     * Dirichlet prior.  This value is not the mode or mean of the Dirichlet but rather the exponential of the expected log weights.
     * note that these effective weights do not add up to 1.  This is fine because in any probabilistic model scaling all weights
     * amounts to an arbitrary normalization constant, but it's important to keep in mind because some classes may expect
     * normalized weights.  In that case the calling code must normalize the weights.
     */
    public double[] effectiveMultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(MathUtils.sum(alpha));
        return Arrays.stream(alpha).map(a -> Math.exp(Gamma.digamma(a) - digammaOfSum)).toArray();
    }

    public int size() { return alpha.length; }
}
