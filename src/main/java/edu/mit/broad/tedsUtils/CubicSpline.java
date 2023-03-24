/*
 * $Id: CubicSpline.java 70611 2008-07-29 18:06:55Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;


/**
 * Cubic splines.  This version only does natural splines on uniform, unit spacing.
 * Adapted from fortran supplied with "A Practical Guide to Splines" by C. de Boor.
 *
 * @author tsharpe
 * @version $Revision: 70611 $
 */
public class CubicSpline
{
	public static void main( String[] args )
	{
	    double[] IDEAL_PEAK = {32., 128., 320., 464., 576., 624., 640., 624., 576., 464., 320., 128., 32.};
	    CubicSpline cs = new CubicSpline(IDEAL_PEAK);
	    System.out.println(cs.argMaxPeak(0,13));
	    System.out.println(cs.argMaxPeak(0,6));
	    System.out.println(cs.argMaxPeak(6,13));
	}

	/**
     * Make a spline interpolator.
     * @param yVals The y values (at x=0, 1, 2...)
     */
    public CubicSpline( double[] yVals )
    {
        this(new DoubleArrayList(yVals));
    }

    /**
     * Make a spline interpolator.
     * @param yVals The y values (at x=0, 1, 2...)
     */
    public CubicSpline( NumberList yVals )
    {
        if ( yVals == null || yVals.size() == 0 )
        {
            throw new IllegalArgumentException("No yVals");
        }
        int len = yVals.size();
        if ( len < 2 )
        {
            throw new IllegalArgumentException("You have to be a consultant to draw a trend line through a single point.");
        }

        mC0 = new double[len];
        for ( int iii = 0; iii < len; ++iii )
        {
            mC0[iii] = yVals.dblVal(iii);
        }
        init(len);
    }

    /**
     * Return the maximum x value supported by the interpolation.
     * (I.e., the number of nodes.)
     */
    public double xMax()
    {
        return mC0.length;
    }

    /**
     * Evaluate the spline interpolation.
     * @param x The x value at which to evaluate the spline.
     * @return The value of the interpolated function.
     */
    public double eval( double x )
    {
        int iii = (int)x;
        if ( iii < 0 )
        {
            iii = 0;
        }
        else if ( iii >= mC0.length )
        {
            iii = mC0.length - 1;
        }
        x -= iii;
        return mC0[iii]+x*(mC1[iii]+x*(mC2[iii]+x*mC3[iii]));
    }

    /**
     * The mean value in a specified range.
     */
    public double meanVal( int start, int end )
    {
        double total = 0.;
        for ( int iii = start; iii < end; ++iii )
        {
            double d = mC0[iii];
            double c = mC1[iii];
            double b = mC2[iii];
            double a = mC3[iii];
            total += a/3 + b/3 + c/2 + d;
        }

        total /= (end - start);
        return total;
    }

    /**
     * Return the x value that will eval to the largest peak in the y values in the specified interval.
     * Note that this isn't the same as the largest y value:  it's the largest y value where the slope = 0.
     * The y value at the start or end of the interval might be larger.
     * If there is no peak in the range, Double.NaN is returned.
     * If there are two or more peaks of equal height, the first is returned.
     * The range searched is inclusive of both the starting point and the ending point.
     * @param start of interval
     * @param end of interval
     * @return x for highest peak y
     */
    public double argMaxPeak( int start, int end )
    {
        if ( start < 0 ) start = 0;
        if ( end >= mC0.length ) end = mC0.length - 1;

        double result = Double.NaN;
        double max = Double.NEGATIVE_INFINITY;
        for ( int iii = start; iii <= end; ++iii )
        {
            double d = mC0[iii];
            double c = mC1[iii];
            double b = mC2[iii];
            double a = mC3[iii];
            double discr = b*b - 3*a*c;
            if ( discr > 0. )
            {
                double minusBOver3a = -b/3/a;
                double discrOver3a = Math.sqrt(discr)/3/a;
                double x = minusBOver3a - discrOver3a;
                if ( x >= 0. && x < 1. && 6*a*x + 2*b < 0 )
                {
                    double y = d+x*(c+x*(b+x*a));
                    if ( y > max )
                    {
                        max = y;
                        result = x + iii;
                    }
                }
                x = minusBOver3a + discrOver3a;
                if ( x >= 0. && x < 1. && 6*a*x + 2*b < 0 )
                {
                    double y = d+x*(c+x*(b+x*a));
                    if ( y > max )
                    {
                        max = y;
                        result = x + iii;
                    }
                }
            }
        }
        return result;
    }

    /**
     * The x value for the largest y value in the range.
     * @param start of range.
     * @param end of range.
     * @return X value for max.
     */
    public double argMax( int start, int end )
    {
        double xVal = argMaxPeak(start,end);
        double yVal = -Double.MAX_VALUE;
        if ( !Double.isNaN(xVal) )
        {
            yVal = eval(xVal);
        }
        double yVal2 = eval(start);
        if ( yVal < yVal2 )
        {
            xVal = start;
            yVal = yVal2;
        }
        yVal2 = eval(end);
        if ( yVal < yVal2 )
        {
            xVal = end;
        }
        return xVal;
    }

    private void init( int len )
    {
        int last = len - 1;

        mC1 = new double[len];
        mC2 = new double[len];
        mC3 = new double[len];

        mC3[0] = 2.;
        for ( int iii = 1; iii < len; ++iii )
        {
            mC3[iii] = (mC0[iii] - mC0[iii - 1]);
        }

        mC1[0] = 3. * mC3[1];
        for ( int iii = 1; iii < last; ++iii )
        {
            double tmp = -1. / mC3[iii-1];
            mC1[iii] = tmp * mC1[iii-1] + 3. * (mC3[iii+1] + mC3[iii]);
            mC3[iii] = tmp + 4.;
        }
        mC1[last] = 3. * mC3[last];
        mC3[last] = 2.;

        double g = -1. / mC3[last - 1];
        mC3[last] += g;
        mC1[last] = (g * mC1[last - 1] + mC1[last]) / mC3[last];

        for ( int iii = last - 1; iii >= 0; --iii )
        {
            mC1[iii] = (mC1[iii] - mC1[iii+1]) / mC3[iii];
        }

        for ( int iii = 0; iii < last; ++iii )
        {
            double divdf1 = (mC0[iii+1] - mC0[iii]);
            double divdf3 = mC1[iii] + mC1[iii+1] - 2. * divdf1;
            mC2[iii] = divdf1 - mC1[iii] - divdf3;
            mC3[iii] = divdf3;
        }
    }

    private double[] mC0;
    private double[] mC1;
    private double[] mC2;
    private double[] mC3;
}
