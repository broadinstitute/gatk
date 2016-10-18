package org.broadinstitute.hellbender.tools.exome.allelicbalancecaller;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;


public class AllelicBalanceCallerModelState implements Serializable {
    private AllelicBalanceCallerModelState() {}

    static final long serialVersionUID = 33733733712L;

    private static final int MIN_CN = 0;
    private static final int MAX_CN = 5;

    private double effectiveAlpha;

    /** Number of rho values we will use.  This implicitly dictates the maximum number of unique values we expect rho
     * to be able to take.
     */
    private int K;

    /** Current unique values of rho in this sample.  For each segment, rho is CCF*purity.  Note that each
     * segment will use a rho from this list.
     */
    private double[] rhos;

    /** The possible absolute copy numbers for allele1.  TYpically, MIN_CN:1:MAX_CN*/
    private int[] mVals;

    /** The possible absolute copy numbers for allele2.  TYpically, MIN_CN:1:MAX_CN*/
    private int[] nVals;

    /** likelihood of CN 0... for all segments.  In other words, all segments share the pis. The entries correspond to
     * the mVals and nVals and will sum to 1.
     */
    private double[] effectivePis;

    /** ploidy * CCF * purity
     *
     * Please note that this does not change in the model.  This is a placeholder.
     */
    private double lambda;

    /** likelihood of each rho.  Will always be of length K */
    private double[] effectivePhis;

    /** The modeled segments, as produced by ACNV */
    private List<ACNVModeledSegment> segments;

    public double[] getEffectivePhis() {
        return effectivePhis;
    }

    public void setEffectivePhis(double[] effectivePhis) {
        this.effectivePhis = effectivePhis;
    }

    public double getEffectiveAlpha() {
        return effectiveAlpha;
    }

    public void setEffectiveAlpha(double eAlpha) {
        this.effectiveAlpha = eAlpha;
    }

    public int getK() {
        return K;
    }

    public void setK(int k) {
        this.K = k;
    }

    public double[] getRhos() {
        return rhos;
    }

    public void setRhos(double[] rhos) {
        this.rhos = rhos;
    }

    public double getLambda() {
        return lambda;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
    }

    public List<ACNVModeledSegment> getSegments() {
        return segments;
    }

    public void setSegments(List<ACNVModeledSegment> segments) {
        this.segments = segments;
    }

    public int[] getmVals() {
        return mVals;
    }

    public void setmVals(int[] mVals) {
        this.mVals = mVals;
    }

    public int[] getnVals() {
        return nVals;
    }

    public void setnVals(int[] nVals) {
        this.nVals = nVals;
    }

    public double[] getEffectivePis() {
        return effectivePis;
    }

    public void setEffectivePis(double[] effectivePis) {
        this.effectivePis = effectivePis;
    }

    /**
     * Create initial state that can be fed into the optimization process.
     * @param rhoThreshold -- minimum rho value to consider.  Anything less than or equal to this number will be disregarded.  Except for zero, that is always considered.
     *                     Must be (0,1)
     * @param segments -- Not {@code null}
     * @param normalPloidy ploidy for a normal sample
     * @return Never {@code null}
     */
    public static AllelicBalanceCallerModelState createInitialCNLOHCallerModelState(final double rhoThreshold, final List<ACNVModeledSegment> segments,
                                                                           final double normalPloidy, final int numRhos) {
        ParamUtils.inRange(rhoThreshold, Math.ulp(0.0), 1.0 - Double.MIN_VALUE, "rho must be (0.0, 1.0)");
        Utils.nonNull(segments);
        ParamUtils.isPositive(numRhos, "Must have at least one rho.");

        AllelicBalanceCallerModelState result = new AllelicBalanceCallerModelState();
        result.setSegments(segments);

        // initialize variables
        result.setEffectiveAlpha(1);
        result.setLambda(normalPloidy);
        result.setK(numRhos);

        result.setEffectivePhis(new double[result.getK()]);
        Arrays.fill(result.getEffectivePhis(), 1.0/result.getK());

        // Initial values of rho.  As a reminder, this will be 0 then skip to the rho threshold and then continue to 1.
        //  There is no point in having initial values of rho that rho is prohibited from taking.
        double[] rhos = new double[]{0};
        rhos = ArrayUtils.addAll(rhos, GATKProtectedMathUtils.createEvenlySpacedPoints(rhoThreshold, 1, result.getK() - 1));
        result.setRhos(rhos);

        result.setmVals(Arrays.stream(GATKProtectedMathUtils.createEvenlySpacedPoints(MIN_CN, MAX_CN, MAX_CN - MIN_CN + 1))
                .mapToInt(d -> (int) Math.round(d)).toArray());
        result.setnVals(Arrays.stream(GATKProtectedMathUtils.createEvenlySpacedPoints(MIN_CN, MAX_CN, MAX_CN - MIN_CN + 1))
                .mapToInt(d -> (int) Math.round(d)).toArray());


        // Change pi mapping in-place.  Set initial values to {0:.01, 1:1, 2:.01, 3:1e-4, ...}
        //  Note this assumes that nVals and mVals are the same.
        result.setEffectivePis(Arrays.stream(result.getmVals()).mapToDouble(m -> Math.pow(1e-2, Math.abs(m-1)))
                .toArray());

        // Put extra emphasis on copy number of 1
        result.getEffectivePis()[1] = 2.0;

        // then normalize to 1, since this is a discrete probability distribution, after all.
        result.setEffectivePis(MathUtils.normalizeFromRealSpace(result.getEffectivePis()));

        return result;
    }
}
