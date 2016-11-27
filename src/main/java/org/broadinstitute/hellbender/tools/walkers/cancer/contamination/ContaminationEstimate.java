package org.broadinstitute.hellbender.tools.walkers.cancer.contamination;

import htsjdk.samtools.util.Locatable;

import java.util.Arrays;

/**
 * a class that estimates and stores the contamination values for a site.
 */
final class ContaminationEstimate {
    private final double precision;             // to what precision do we want to run; e.g. if set to 1, we run using 1% increments
    private final double[] bins;                // the bins representing the discrete contamination levels we're evaluating
    private double populationFit = 0.0;
    private String popultationName = "";

    private static double[] precalculatedEpsilon;

    private int arrayAlleleObservations = 0;
    private int alternateAlleleObservations = 0;

    // precalculate the 128 values of epsilon that are possible
    static {
        precalculatedEpsilon = new double[Byte.MAX_VALUE+1];

        for(int i=0; i <= (int)Byte.MAX_VALUE; i++) {
            precalculatedEpsilon[i] = Math.pow(10.0,-1.0*(((double)i)/10.0));
        }
    }

    /**
     * create the contamination estimate, given:
     * @param precision the precision value, to what level are we calculating the contamination
     */
    public ContaminationEstimate(double precision,
                                 double maf,
                                 byte[] bases,
                                 byte[] quals,
                                 byte arrayAllele,
                                 byte hapmapAlt,
                                 String popName,
                                 Locatable locus
    ) {
        // setup the bins to the correct precision
        this.precision = precision;
        bins = new double[(int)Math.ceil(100/precision)+1];
        if (maf == 0) maf = 0.00001;

        popultationName = popName;

        Arrays.fill(bins,0.0); // just to make sure we don't have any residual values

        // convert the quals
        double[] realQuals = new double[quals.length];
        int qIndex = 0;
        for (byte qual : quals) {realQuals[qIndex++] = Math.pow(10.0,-1.0*(qual/10.0));}

        // check our inputs
        if (maf > 1.0 || maf < 0.0) throw new IllegalArgumentException("Invalid allele Freq: must be between 0 and 1 (inclusive), maf was " + maf + " for population " + popName);

        // calculate the contamination for each bin
        int qualOffset = 0;
        for (byte base : bases) {

            if (base == arrayAllele) { arrayAlleleObservations++; }
            if (base == hapmapAlt) { alternateAlleleObservations++; }
            double epsilon = precalculatedEpsilon[quals[qualOffset++]];

            for (int index = 0; index < bins.length; index++) {


                double contaminationRate = (1.0 - (double) index / (double) bins.length);

                if (base == arrayAllele) {
                    bins[index] += Math.log((1.0 - contaminationRate) * (1.0 - epsilon) +
                            contaminationRate * ((maf) * (1.0 - epsilon) + (1.0 - maf) * (epsilon/3.0)));
                    populationFit += Math.log(epsilon);

                } else if(hapmapAlt == base) {
                    bins[index] += Math.log((1.0 - contaminationRate) * (epsilon / 3.0) +
                            contaminationRate * ((maf) * (epsilon/3.0) + (1.0 - maf) * (1.0 - epsilon)));

                    populationFit += Math.log(maf + epsilon);
                }
            }
        }
    }

    public double[] getBins() {
        return bins;
    }

    public void setPopulationFit(double populationFit) {
        this.populationFit = populationFit;
    }

    public double getPopulationFit() {
        return populationFit;
    }

    public String getPopultationName() {
        return popultationName;
    }

    public static class ConfidenceInterval {

        private double start;
        private double stop;
        private double contamination;
        private double maxLikelihood;
        double[] newBins;

        public ConfidenceInterval(double bins[], double intervalArea) {
            // make a copy of the bins in non-log space
            int maxIndex = 0;
            for (int x = 0; x < bins.length; x++) if (bins[x] > bins[maxIndex]) maxIndex = x;
            newBins = new double[bins.length];
            maxLikelihood = bins[maxIndex];

            int index = 0;
            double total = 0.0;
            for (double d : bins) {
                newBins[index] = Math.pow(10,(bins[index] - bins[maxIndex]));
                total += newBins[index];
                index++;
            }

            for (int x = 0; x < newBins.length; x++) {
                newBins[x] = newBins[x] / total;
            }
            double areaUnderCurve = 0;
            int leftIndex = maxIndex;
            int rightIndex = maxIndex;
            while (areaUnderCurve < 0.95) {

                // if the "left" bin is bigger, and can be moved, move it
                if (newBins[leftIndex] >= newBins[rightIndex] && leftIndex > 0) {
                    leftIndex--;
                } else {
                    // otherwise move the right bin if possible
                    if (rightIndex < bins.length - 1) {
                        rightIndex++;
                    } else {
                        // and if not move the left bin, or die
                        if (leftIndex > 0) {
                            leftIndex--;
                        } else {
                            throw new RuntimeException("Error trying to compute confidence interval");
                        }
                    }
                }

                areaUnderCurve = 0.0;
                for (int x = leftIndex; x <= rightIndex; x++)
                    areaUnderCurve += newBins[x];
            }
            start = (bins.length - rightIndex) * (100.0/bins.length);
            stop = (bins.length - leftIndex) * (100.0/bins.length);
            contamination = (bins.length - maxIndex) * (100.0/bins.length);
        }

        public double getStart() {
            return start;
        }

        public double getStop() {
            return stop;
        }

        public double getContamination() {
            return contamination;
        }

        public double getMaxLikelihood() {
            return maxLikelihood;
        }

        public String toString() {
            return contamination + "[" + start + " - " + stop + "] log likelihood = " + maxLikelihood;
        }
    }
}
