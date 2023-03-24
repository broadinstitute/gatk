/*
 * $Id: MathFunctions.java 74943 2008-11-07 22:10:54Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

/**
 * Some handy statistics functions. Most of the stuff in here comes from:
 * 
 * <pre>
 * 
 *  Numerical Recipes in C
 *  The Art of Scientific Computing
 *  by Press, Flannery, Teukolsky, and Vetterling
 *  (c) Cambridge University Press 1988
 *  
 * </pre>
 * 
 * @author Ted Sharpe
 * @version $Revision: 16386 $
 */
public class MathFunctions
{
    public static final double EPS = calcEpsilon();
    public static final double LN_EPS = Math.log(calcEpsilon());

    /**
     * Fisher's z transformation.
     */
    public static double fisherZ( double r )
    {
        return Math.log((1. + r) / (1. - r)) / 2;
    }

    /**
     * Caclulate partial correlation.
     */
    public static double partialCorrelation( double r01, double r0c, double r1c )
    {
        return (r01 - r0c * r1c) / (Math.sqrt(1. - r0c * r0c) * Math.sqrt(1. - r1c * r1c));
    }

    /**
     * Calculates the cumulative distribution function for Spearman's rank correlation.
     */
    public static double rsCDF( double rs, int nnn )
    {
        return tCDF(tStat(rs, nnn), nnn - 2);
    }

    /**
     * Calculates the t value for rank correlation.
     */
    public static double tStat( double rs, int nnn )
    {
        return rs * Math.sqrt((nnn - 2) / (1. - rs * rs));
    }

    /**
     * Calculates the cumulative distribution function of Student's t distrubtion.
     */
    public static double tCDF( double t, int df )
    {
        return MathFunctions.betai(df / 2., .5, df / (df + t * t));
    }

    /**
     * Calculates the incomplete beta function.
     */
    public static double betai( double a, double b, double x )
    {
        if ( x < 0. || x > 1. )
            throw new IllegalArgumentException("Out of range: " + x);

        double result = 0.;
        if ( x != 0. && x != 1. )
            result = Math.exp(lnGamma(a + b) - lnGamma(a) - lnGamma(b) + a * Math.log(x) + b * Math.log(1. - x));

        if ( x < (a + 1.) / (a + b + 2.) )
            result *= betacf(a, b, x) / a;
        else
            result = 1. - result * betacf(b, a, 1. - x) / b;

        return result;
    }

    private static double betacf( double a, double b, double x )
    {
        int m, m2;
        double aa, del;

        double qab = a + b;
        double qap = a + 1.;
        double qam = a - 1.;
        double c = 1.;
        double d = 1. - qab * x / qap;
        if ( Math.abs(d) < MIN_VALUE )
            d = MIN_VALUE;
        d = 1. / d;
        double h = d;

        for ( m = 1; m <= MAX_ITR; ++m )
        {
            m2 = 2 * m;
            aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1. + aa * d;
            if ( Math.abs(d) < MIN_VALUE )
                d = MIN_VALUE;
            c = 1. + aa / c;
            if ( Math.abs(c) < MIN_VALUE )
                c = MIN_VALUE;
            d = 1. / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1. + aa * d;
            if ( Math.abs(d) < MIN_VALUE )
                d = MIN_VALUE;
            c = 1. + aa / c;
            if ( Math.abs(c) < MIN_VALUE )
                c = MIN_VALUE;
            d = 1. / d;
            del = d * c;
            h *= del;
            if ( Math.abs(del - 1.) < EPS )
                return h;
        }

        throw new IllegalStateException("failed to converge after " + MAX_ITR + " iterations");
    }

    /**
     * Computes the cumulative distribution function of the gamma distribution.
     * Keeps track of the shape and scale params so that you can evaluate easily.
     * @author tsharpe
     */
    public static final class GammaDist
    {
        public GammaDist( double shape, double scale )
        {
            mShape = shape;
            mScale = scale;
        }

        public double mean()
        {
            return mShape*mScale;
        }

        public double mode()
        {
            return (mShape-1)*mScale;
        }

        public double stdDev()
        {
            return Math.sqrt(mShape)*mScale;
        }

        public double evalCDF( double x )
        {
            return 1. - incompleteGamma(mShape,x/mScale);
        }

        private double mShape;
        private double mScale;
    }

    /**
     * This calculates the incompleteGamma function Q(a,x) = 1 - P(a,x) = littleGamma(a,x)/Gamma(a).
     */
    public static double incompleteGamma( double a, double x )
    {
        double result = Double.NaN;

        if ( x >= 0. && a > 0. )
        {
            double tmp = Math.exp(a * Math.log(x) - x - lnGamma(a));

            if ( x < a + 1. )
            {
                double del = 1. / a;
                double sum = del;

                do
                {
                    a += 1.;
                    del *= x / a;
                    sum += del;
                }
                while ( del > sum * EPS );

                result = 1. - sum * tmp;
            }
            else
            {
                double a0 = 1.;
                double a1 = x;
                double an = 0.;
                double b0 = 0.;
                double b1 = 1.;
                double fac = 1.;
                double lastG = 0.;

                while ( true )
                {
                    an += 1.;
                    double ana = an - a;
                    a0 = (a1 + a0 * ana) * fac;
                    b0 = (b1 + b0 * ana) * fac;
                    double anf = an * fac;
                    a1 = x * a0 + anf * a1;
                    b1 = x * b0 + anf * b1;
                    if ( a1 != 0. )
                    {
                        fac = 1. / a1;
                        double g = b1 * fac;
                        if ( Math.abs((g - lastG) / g) < EPS )
                        {
                            result = tmp * g;
                            break;
                        }
                        lastG = g;
                    }
                }
            }
        }

        return result;
    }

    /**
     * This calculates the natural logarithm of the gamma function. (The one with a single argument and spelled with the
     * capital version of the third letter of the Greek alphabet that looks like what you draw at the beginning of a
     * game of hangman.) The method is based on an approximation derived by Lanczos. This approximation works well for
     * arguments > 1. We do a fix-up for arguments > 0., and though the gamma function is defined for all reals except
     * non-positive integers, we're just going to return NaN for arguments <= 0, because we don't need them for
     * calculating the cumulative distribution function of chi squared, which is why I typed this in.
     * 
     * @param x
     *            x
     * @return ln(Gamma(x))
     */
    public static double lnGamma( double x )
    {
        if ( x <= 0. )
            throw new IllegalArgumentException("out of range");

        return -5.5 - x + (x + .5) * Math.log(x + 5.5) +
                Math.log(2.5066282746310005 / x * (1.000000000190015 + 76.18009172947146 / (x + 1.) - 86.50532032941677 / (x + 2.)
                             + 24.01409824083091 / (x + 3.) - 1.231739572450155 / (x + 4.) + .1208650973866179e-2 / (x + 5.) - .5395239384953e-5 / (x + 6.)));
    }

    /**
     * Calculates the complementary error function. This Chebyshev approximation is good to just shy of 7 digits.
     */
    public static double erfc( double x )
    {
        double z = Math.abs(x);
        double t = 1. / (1. + z / 2.);
        double result = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
        return x >= 0. ? result : 2. - result;
    }

    /**
     * Calculates the cumulative distribution function of the unit normal distribution.
     * I.e., it's the probability that a z score is less than the argument.
     * I.e., it's the integral of the unit normal distribution from -infinity to z.
     * Uses the erfc implementation above, so it's only good to 6 or 7 digits.
     */
    public static double zCDF( double z )
    {
        return (2. - erfc(z / SQRT_2)) / 2.;
    }

    /**
     * Calculates the likelihood of a z score lesser or equal in magnitude than the one presented.
     * I.e., it's the absolute value of the integral of the unit normal distribution from -z to z.
     * Uses the erfc implementation above, so it's only good to 6 or 7 digits.
     */
    public static double zOrLessProbability( double z )
    {
        return 1. - erfc(Math.abs(z) / SQRT_2);
    }

    /**
     * Calculates the likelihood of a z score greater or equal in magnitude than the one presented.
     * I.e., it's 1. less the absolute value of the integral of the unit normal distribution from -z to z.
     * Uses the erfc implementation above, so it's only good to 6 or 7 digits.
     */
    public static double zOrMoreProbability( double z )
    {
        return erfc(Math.abs(z) / SQRT_2);
    }

    /**
     * Calculate log(a+b) given log(a) and log(b)
     */
    public static double sumLogs( double logA, double logB )
    {
        double result;

        if ( logA < logB )
            result = sumLogs(logB, logA);
        else
        {
            double diff = logB - logA;
            if ( diff < LN_EPS )
                result = logA;
            else
                result = logA + Math.log(1. + Math.exp(diff));
        }

        return result;
    }

    /**
     * Calculates the machine epsilon.
     */
    private static double calcEpsilon()
    {
        double eps = 1.;

        while ( 1. + eps > 1. )
            eps /= 2.;

        return eps;
    }

    private static final int MAX_ITR = 100;
    private static final double MIN_VALUE = 10. * Double.MIN_VALUE;
    private static final double SQRT_2 = Math.sqrt(2.);
}
