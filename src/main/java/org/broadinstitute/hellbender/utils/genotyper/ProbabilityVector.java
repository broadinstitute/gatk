package org.broadinstitute.hellbender.utils.genotyper;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;

public class ProbabilityVector {
    private final double[] probabilityArray;
    private final int minVal;
    private final int maxVal;

    final static double LOG_DYNAMIC_RANGE = 10; // values X below max vector value will be removed

    /**
     * Default constructor: take vector in log-space, with support from range [0,len-1]
     * @param vec                                  Probability (or likelihood) vector in log space
     * @param compressRange                        If true, compress by eliminating edges with little support
     */
    public ProbabilityVector(final double[] vec, final boolean compressRange) {

        final int maxValIdx = MathUtils.maxElementIndex(vec);
        final double maxv = vec[maxValIdx];
        if (maxv > 0.0) {
            throw new GATKException("BUG: Attempting to create a log-probability vector with positive elements");
        }

        if (compressRange) {
            minVal = getMinIdx(vec, maxValIdx);
            maxVal = getMaxIdx(vec, maxValIdx);
            probabilityArray = Arrays.copyOfRange(vec, minVal, maxVal + 1);

        }   else {
            probabilityArray = vec;
            minVal = 0;
            maxVal = vec.length-1;

        }
    }

    public ProbabilityVector(final double[] vec) {
        this(vec,true);
    }

    public ProbabilityVector(final ProbabilityVector other, final boolean compressRange) {
        // create new probability vector from other.
        this(other.getUncompressedProbabilityVector(), compressRange);
        
    }
    public int getMinVal() { return minVal;}
    public int getMaxVal() { return maxVal;}
    public double[] getProbabilityVector() { return probabilityArray;}
    
    public double[] getProbabilityVector(final int minVal, final int maxVal) {
        // get vector in specified range. If range is outside of current vector, fill with negative infinities
        final double[] x = new double[maxVal - minVal + 1];

        for (int k=minVal; k <= maxVal; k++) {
            x[k - minVal] = getLogProbabilityForIndex(k);
        }


        return x;
    }

    public double[] getUncompressedProbabilityVector() {
        final double[] x = new double[maxVal+1];
        
        for (int i=0; i < minVal; i++) {
            x[i] = Double.NEGATIVE_INFINITY;
        }
        for (int i=minVal; i <=maxVal; i++) {
            x[i] = probabilityArray[i - minVal];
        }

        return x;
    }
    /**
     * Return log Probability for original index i
     * @param idx   Index to probe
     * @return      log10(Pr X = i) )
     */
    public double getLogProbabilityForIndex(final int idx) {
    if (idx < minVal || idx > maxVal) {
        return Double.NEGATIVE_INFINITY;
    } else {
        return probabilityArray[idx - minVal];
    }
    }

    //public ProbabilityVector
    public static ProbabilityVector compressVector(final double[] vec ) {
        return new ProbabilityVector(vec, true);
    }

    /**
     * Determine left-most index where a vector exceeds (max Value - DELTA)
     * @param vec                    Input vector
     * @param maxValIdx              Index to stop - usually index with max value in vector
     * @return                       Min index where vector > vec[maxValIdx]-LOG_DYNAMIC_RANGE
     */
    private static int getMinIdx(final double[] vec, final int maxValIdx) {
        int edgeIdx;
        for (edgeIdx=0; edgeIdx<=maxValIdx; edgeIdx++ ) {
            if (vec[edgeIdx] > vec[maxValIdx]-LOG_DYNAMIC_RANGE) {
                break;
            }
        }

        return edgeIdx;


    }

    /**
     * Determine right-most index where a vector exceeds (max Value - DELTA)
     * @param vec                    Input vector
     * @param maxValIdx              Index to stop - usually index with max value in vector
     * @return                       Max index where vector > vec[maxValIdx]-LOG_DYNAMIC_RANGE
     */
    private static int getMaxIdx(final double[] vec, final int maxValIdx) {
        int edgeIdx;
        for (edgeIdx=vec.length-1; edgeIdx>=maxValIdx; edgeIdx-- ) {
            if (vec[edgeIdx] > vec[maxValIdx]-LOG_DYNAMIC_RANGE) {
                break;
            }
        }

        return edgeIdx;


    }

    /**
     *
     * @param other
     * @return
     */
    public double logDotProduct(final ProbabilityVector other) {
        // find overlap in range
        final int minRange = Math.max(this.minVal, other.getMinVal());
        final int maxRange = Math.min(this.maxVal, other.getMaxVal());
        if (minRange > maxRange) {
            return Double.NEGATIVE_INFINITY;
        }

        // x = 0,1,2,   y = 2,3,4. minRange = 2, maxRange = 2
        final double[] result = new double[maxRange - minRange+1];
        for (int k=0; k <= maxRange-minRange; k++) {
            final int startI = minRange - this.minVal;
            final int startJ = minRange - other.getMinVal();
            result[k] = this.probabilityArray[k+startI] + other.probabilityArray[k+startJ];
            


        }
        return MathUtils.approximateLog10SumLog10(result);
    }

}
