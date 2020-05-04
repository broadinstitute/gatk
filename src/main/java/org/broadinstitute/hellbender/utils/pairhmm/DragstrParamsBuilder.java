package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.commons.lang.math.IntRange;
import org.broadinstitute.hellbender.utils.Utils;

public class DragstrParamsBuilder {

    private final int maxPeriod;
    private final int maxRepeats;

    public final double[][] gp;
    public final double[][] gcp;
    public final double[][] api;

    public DragstrParamsBuilder(final int maxPeriod, final int maxRepeats) {
        this.maxPeriod = maxPeriod;
        this.maxRepeats = maxRepeats;
        this.gp =  new double[maxPeriod][maxRepeats];
        this.gcp = new double[maxPeriod][maxRepeats];
        this.api = new double[maxPeriod][maxRepeats];
    }

    public void set(final int period, final IntRange repeatRange, final double gp, final double gcp, final double api) {
        for (int i = repeatRange.getMinimumInteger(); i <= repeatRange.getMaximumInteger(); i++) {
            set(period, i, gp, gcp, api);
        }
    }

    public void set(final int period, final int repeats, final double gp, final double gcp, final double api) {
        final int p = period - 1;
        final int r = repeats - 1;
        Utils.validIndex(p, maxPeriod);
        Utils.validIndex(r, maxRepeats);
        this.gp[p][r] = gp;
        this.gcp[p][r] = gcp;
        this.api[p][r] = api;
    }

    public DragstrParams make(final double minGop, final double maxGop, final double stepwise) {
        final double[][] gop = new double[gp.length][gp[0].length];
        for (int i = 0; i < gop.length; i++) {
            for (int j = 0; j < gop[i].length; j++) {
                gop[i][j] = gopCalculation(gp[i][j], gcp[i][j], i, minGop, maxGop, stepwise);
            }
        }
        return DragstrParams.of(maxPeriod, maxRepeats, gop, gcp, api);
    }

    private double gopCalculation(final double gp, final double gcp, final int period, final double minGop, final double maxGop, final double stepwise) {

        double gpProb = Math.pow(10, -.1 * gp);
        double bestGop = minGop;
        final double c = Math.pow(10, -0.1 * gcp);
        double bestCost = Double.POSITIVE_INFINITY;
        for (double gop = minGop; maxGop - gop > -0.001; gop += stepwise) {
            final double g = Math.pow(10, -0.1 * gop);
            final double p_gap = g * Math.pow(c, period - 1) * (1.0 - c);
            final double p_no_gap = Math.pow(1 - 2 * g, period + 1);
            final double cost = Math.abs(p_gap / p_no_gap - gpProb);
            if (cost < bestCost) {
                bestGop = gop;
                bestCost = cost;
            }
        }
        return bestGop;
    }
}
