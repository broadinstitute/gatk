package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;

/**
 * Represents a combination of STR model parameters pi and lambda together with the training size.
 */
public class STRModelParameters {

    private static final double LN_10 = Math.log(10);
    public double modelFreeLikelihood = 0;
    public double likelihood = 0;
    public int iterations = 0;
    public final double pi;
    public final double lambda;
    private final double expOfMinusLambda;
    private final double oneMinusExpOfMinusLambda;
    private final double inverseLambda;
    private final double inverseOneMinusExpOfMinusLambda;
    public final long size; // number of cases used to calculate the lambda and error-rate.
    private final double log10OffDiagonalFixedFactor;
    private final double log10Pi;
    private final double logLambda;
    private final double log10OneMinusPi;

    STRModelParameters(final long size, final double pi, final double lambda) {
        this.size = size;
        this.lambda = lambda;
        this.logLambda = Math.log(lambda);
        this.expOfMinusLambda = Math.exp(-lambda);
        this.inverseLambda = 1.0 / lambda;
        this.oneMinusExpOfMinusLambda = 1 - expOfMinusLambda;
        this.inverseOneMinusExpOfMinusLambda = 1.0 / oneMinusExpOfMinusLambda;
        this.log10Pi = Math.log10(pi);
        this.log10OffDiagonalFixedFactor = log10Pi - MathUtils.log10(2) - MathUtils.log10OneMinusX(expOfMinusLambda);
        this.log10OneMinusPi = MathUtils.log10OneMinusX(pi);
        this.pi = pi;
    }

    public static STRModelParameters parse(final String spec) {
        final String[] parts = spec.split("\\s+");
        if (parts.length < 3) {
            throw new IllegalArgumentException("Bad spec, does not have at least d3 values");
        } else {
            final int lastIndex = parts.length - 1;
            return new STRModelParameters(Integer.parseInt(parts[lastIndex - 2]), Double.parseDouble(parts[lastIndex - 1]), Double.parseDouble(parts[lastIndex]));
        }
    }

    RealMatrix log10TransformationMatrix(final List<STRAllele> alleles) {
        final Array2DRowRealMatrix result = new Array2DRowRealMatrix(alleles.size(), alleles.size());
        final int alleleCount = alleles.size();

        for (int i = 0; i < alleleCount; i++) {
            result.setEntry(i, i, log10OneMinusPi);
            for (int j = i + 1; j < alleleCount; j++) {
                final int rcd = Math.abs(alleles.get(i).repeatCount - alleles.get(j).repeatCount);
                final double value = (rcd * logLambda - lambda - MathUtils.log10Factorial(rcd)) + log10OffDiagonalFixedFactor;
                result.setEntry(i, j, value);
                result.setEntry(j, i, value);
            }
        }
        
        return result;
    }

    /**
     * Notice that the result derivative matrix is not log10ed.
     * @param alleles
     * @return never {@code null}.
     */
    RealMatrix lambdaDerivativeMatrix(final List<STRAllele> alleles) {
        final int alleleCount = alleles.size();
        final Array2DRowRealMatrix result = new Array2DRowRealMatrix(alleleCount, alleleCount);
        final double offDiagonalSecondPartValue = oneMinusExpOfMinusLambda;
        final double offDiagonalSecondPartDerivate = expOfMinusLambda;
        final double offDiagonalDenominatorInverse = 1.0 / (offDiagonalSecondPartValue * offDiagonalSecondPartValue);
        for (int i = 0; i < alleleCount; i++) {
            result.setEntry(i, i, 0);
            for (int j = i + 1; j < alleleCount; j++) {
                final int rcd = Math.abs(alleles.get(i).repeatCount - alleles.get(j).repeatCount);
                final double scale = .5 * pi * Math.exp(-MathUtils.log10Factorial(rcd) * LN_10);
                final double lambdaPowRcd = Math.pow(lambda, rcd);
                final double offDiagonalFirstPartValue = lambdaPowRcd * expOfMinusLambda;
                final double offDiagonalFirstPartDerivate = (rcd * lambdaPowRcd * inverseLambda - lambdaPowRcd) * expOfMinusLambda;
                final double value = scale * ((offDiagonalFirstPartDerivate * offDiagonalSecondPartValue)
                        - (offDiagonalFirstPartValue * offDiagonalSecondPartDerivate)) * offDiagonalDenominatorInverse;
                result.setEntry(i, j, value);
                result.setEntry(j, i, value);
            }
        }
        return result;
    }

    /**
     * Notice that the result derivative matrix is not log10ed.
     * @param alleles
     * @return never {@code null}.
     */
    RealMatrix piDerivativeMatrix(final List<STRAllele> alleles) {
        final int alleleCount = alleles.size();
        final Array2DRowRealMatrix result = new Array2DRowRealMatrix(alleleCount, alleleCount);
        for (int i = 0; i < alleleCount; i++) {
            result.setEntry(i, i, -1.0);
            for (int j = i + 1; j < alleleCount; j++) {
                final int rcd = Math.abs(alleles.get(i).repeatCount - alleles.get(j).repeatCount);
                final double value = .5 * Math.exp(rcd * logLambda - lambda - MathUtils.log10Factorial(rcd) * LN_10) * inverseOneMinusExpOfMinusLambda;
                result.setEntry(i, j, value);
                result.setEntry(j, i, value);
            }
        }
        return result;
    }

    public String toString() {
        return String.format("[STR-Param-Set size=%d pi=%g lambda=%g iterations=%d likelihood=%g initialLikelihood=%g]", size, pi, lambda, iterations, likelihood, modelFreeLikelihood);
    }

    public VariantContext annotate(final VariantContext original) {
        final VariantContextBuilder builder = new VariantContextBuilder(original);
        builder.attribute(STRModel.PI_INFO_KEY, String.format("%.4g", pi));
        builder.attribute(STRModel.TAU_INFO_KEY, String.format("%.4g", lambda));
        return builder.make();
    }
}
