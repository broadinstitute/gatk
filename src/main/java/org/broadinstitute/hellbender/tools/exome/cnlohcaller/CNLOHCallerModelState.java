package org.broadinstitute.hellbender.tools.exome.cnlohcaller;

import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.IntStream;


public class CNLOHCallerModelState implements Serializable {
    private CNLOHCallerModelState() {}

    static final long serialVersionUID = 33733733712L;

    final static private int MIN_CN = 0;
    final static private int MAX_CN = 5;

    private double eAlpha;

    /** Number of rho values we will use.  This implicitly dictates the maximum number of unique values we expect rho
     * to be able to take.
     */
    private int numK;

    /** Current unique values of rho in this sample.  For each segment, rho is CCF*purity.  Note that each
     * segment will use a rho from this list.
     */
    private double[] rhos;

    /** likelihood of CN 0... for all segments.  In other words, all segments share the pis. The keys are the actual
     * copy numbers we are allowing.  The values are the pi (likelihood of that copy number).
     *
     * A SortedMap is best for this map.
     */
    private SortedMap<Integer, Double> cnToPiMap = new TreeMap<>();

    /** ploidy * CCF * purity
     *
     * Please note that this does not change in the model.  This is a placeholder.
     */
    private double lambda;

    /** Log'ed likelihood of each rho.  Will always be of length K */
    private double[] ELnPhiK;

    /** The modeled segments, as produced by ACNV */
    private List<ACNVModeledSegment> segments;


    public double[] getELnPhiK() {
        return ELnPhiK;
    }

    public void setELnPhiK(double[] ELnPhiK) {
        this.ELnPhiK = ELnPhiK;
    }

    public double getEAlpha() {
        return eAlpha;
    }

    public void setEAlpha(double eAlpha) {
        this.eAlpha = eAlpha;
    }

    public int getNumK() {
        return numK;
    }

    public void setNumK(int numK) {
        this.numK = numK;
    }

    public double[] getRhos() {
        return rhos;
    }

    public void setRhos(double[] rhos) {
        this.rhos = rhos;
    }

    public SortedMap<Integer, Double> getCnToPiMap() {
        return cnToPiMap;
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

    /**
     * Create initial state that can be fed into the optimization process.
     * @param rhoThreshold -- minimum rho value to consider.  Anything less than or equal to this number will be disregarded.  Except for zero, that is always considered.
     *                     Must be (0,1)
     * @param segments -- Not {@code null}
     * @return Never {@code null}
     */
    public static CNLOHCallerModelState createInitialCNLOHCallerModelState(final double rhoThreshold, final List<ACNVModeledSegment> segments) {

        ParamUtils.inRange(rhoThreshold, Math.ulp(0.0), 1.0 - Double.MIN_VALUE, "rho must be (0.0, 1.0)");
        Utils.nonNull(segments);

        CNLOHCallerModelState result = new CNLOHCallerModelState();

        // initialize variables
        result.setEAlpha(1);
        result.setLambda(2.0);
        result.setNumK(25);

        result.setELnPhiK(new double[result.getNumK()]);
        Arrays.fill(result.getELnPhiK(), Math.log(1.0/result.getNumK()));

        // Initial values of rho.  As a reminder, this will be 0 then skip to the rho threshold and then continue to 1.
        //  There is no point in having initial values of rho that rho is prohibited from taking.
        double[] rhos = new double[]{0};
        rhos = ArrayUtils.addAll(rhos, GATKProtectedMathUtils.linspace(rhoThreshold, 1, result.getNumK() - 1));
        result.setRhos(rhos);

        // Change pi mapping in-place.  Set initial values to {0:.01, 1:1, 2:.01, 3:1e-4, ...}
        Map<Integer, Double> cnToPiMap = result.getCnToPiMap();
        IntStream.range(CNLOHCallerModelState.MIN_CN, CNLOHCallerModelState.MAX_CN+1).
                forEach(cn -> cnToPiMap.put(cn, Math.pow(1e-2, Math.abs(cn-1))));
        cnToPiMap.put(1, 2.0);

        // then normalize to 1, since this is a discrete probability distribution, after all.
        final double piTotal = cnToPiMap.values().stream().reduce(0.0, Double::sum);
        cnToPiMap.keySet().stream().forEach(k -> cnToPiMap.put(k, cnToPiMap.get(k)/piTotal));

        result.setSegments(segments);

        return result;
    }
}
