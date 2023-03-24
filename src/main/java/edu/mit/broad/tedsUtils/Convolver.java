/**
 * $Id: Convolver.java 72836 2008-09-19 15:56:29Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils;


/**
 * Calculates a convolution.
 * The size of the input and output is the same.  This is achieved by
 * ignoring the leading K/2 (where K is the length of the kernel) and
 * trailing K/2 outputs.  These outputs aren't terribly useful, because
 * the input signal is less than half enmeshed in the kernel at these
 * points.  An alternative point of view:  Only those samples where the
 * input is exposed to the central point of an odd-length kernel will
 * be output, which is handy for zero-phase filters.
 *
 * @author tsharpe
 * @version $Revision: 72836 $
 */
public class Convolver
{
    public Convolver( double[] kernel )
    {
        mKernel = kernel.clone();
    }

    /**
     * Scales the kernel so that it's sum is equal to the specified gain.
     */
    public void setGain( double gain )
    {
        int len = mKernel.length;
        double tot = 0.;
        for ( int idx = 0; idx < len; ++idx )
        {
            tot += mKernel[idx];
        }
        double scale = gain / tot;
        for ( int idx = 0; idx < len; ++idx )
        {
            mKernel[idx] *= scale;
        }
    }

    /**
     * Perform the convolution.
     */
    public NumberList convolve( NumberList input )
    {
        int kernelLen = mKernel.length;
        int offset = kernelLen / 2;
        int inputLen = input.size();
        double[] outputs = new double[inputLen];
        for ( int iii = 0; iii < inputLen; ++iii )
        {
            int outputIdx = iii - offset;
            int jStop = Math.min(kernelLen,inputLen - outputIdx);
            int jStart = 0;
            if ( outputIdx < 0 )
            {
                jStart = -outputIdx;
                outputIdx = 0;
            }
            double iVal = input.dblVal(iii);
            while ( jStart < jStop )
            {
                outputs[outputIdx++] += iVal * mKernel[jStart++];
            }
        }
        return new DoubleArrayList(outputs);
    }

    private double[] mKernel;
}
