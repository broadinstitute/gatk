package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class TrancheManager {

    protected final static Logger logger = LogManager.getLogger(TrancheManager.class);

    // ---------------------------------------------------------------------------------------------------------
    //
    // Code to determine FDR tranches for VariantDatum[]
    //
    // ---------------------------------------------------------------------------------------------------------

    public static abstract class SelectionMetric {
        String name = null;

        public SelectionMetric(String name) {
            this.name = name;
        }

        public String getName() { return name; }

        public abstract double getThreshold(double tranche);
        public abstract double getTarget();
        public abstract void calculateRunningMetric(List<VariantDatum> data);
        public abstract double getRunningMetric(int i);
        public abstract int datumValue(VariantDatum d);
    }

    public static class TruthSensitivityMetric extends SelectionMetric {
        double[] runningSensitivity;
        int nTrueSites = 0;

        public TruthSensitivityMetric(final int nTrueSites) {
            super("TruthSensitivity");
            this.nTrueSites = nTrueSites;
        }

        public double getThreshold(final double tranche) {
            return 1.0 - tranche/100.0; // tranche of 1 => 99% sensitivity target
        }

        public double getTarget() { return 1.0; }

        public void calculateRunningMetric(final List<VariantDatum> data) {
            int nCalledAtTruth = 0;
            runningSensitivity = new double[data.size()];

            for ( int i = data.size() - 1; i >= 0; i-- ) {
                VariantDatum datum = data.get(i);
                nCalledAtTruth += datum.atTruthSite ? 1 : 0;
                runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
            }
        }

        public double getRunningMetric(final int i) {
            return runningSensitivity[i];
        }

        public int datumValue(final VariantDatum d) {
            return d.atTruthSite ? 1 : 0;
        }
    }

    public static List<Tranche> findTranches( final List<VariantDatum> data,
                                              final List<Double> tranches,
                                              final SelectionMetric metric,
                                              final VariantRecalibratorArgumentCollection.Mode model ) {
        return findTranches( data, tranches, metric, model, null );
    }

    public static List<Tranche> findTranches(
            final List<VariantDatum> data,
            final List<Double> trancheThresholds,
            final SelectionMetric metric,
            final VariantRecalibratorArgumentCollection.Mode model,
            final File debugFile ) {
        logger.info(String.format("Finding %d tranches for %d variants", trancheThresholds.size(), data.size()));

        Collections.sort( data, VariantDatum.VariantDatumLODComparator );
        metric.calculateRunningMetric(data);

        if ( debugFile != null) {
            writeTranchesDebuggingInfo(debugFile, data, metric);
        }

        List<Tranche> tranches = new ArrayList<>();
        for ( double trancheThreshold : trancheThresholds ) {
            Tranche t = findTranche(data, metric, trancheThreshold, model);

            if ( t == null ) {
                if ( tranches.size() == 0 ) {
                    throw new UserException(String.format(
                            "Couldn't find any tranche containing variants with a %s > %.2f. Are you sure the truth files contain unfiltered variants which overlap the input data?",
                            metric.getName(),
                            metric.getThreshold(trancheThreshold)));
                }
                break;
            }

            tranches.add(t);
        }

        return tranches;
    }

    private static void writeTranchesDebuggingInfo(final File f, final List<VariantDatum> tranchesData, final SelectionMetric metric ) {
        try {
            PrintStream out = new PrintStream(f);
            out.println("Qual metricValue runningValue");
            for ( int i = 0; i < tranchesData.size(); i++ ) {
                VariantDatum  d = tranchesData.get(i);
                int score = metric.datumValue(d);
                double runningValue = metric.getRunningMetric(i);
                out.printf("%.4f %d %.4f%n", d.lod, score, runningValue);
            }
            out.close();
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(f, e);
        }
    }

    public static Tranche findTranche(
            final List<VariantDatum> data,
            final SelectionMetric metric,
            final double trancheThreshold,
            final VariantRecalibratorArgumentCollection.Mode model ) {
        logger.info(String.format("  Tranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, metric.getThreshold(trancheThreshold)));

        double metricThreshold = metric.getThreshold(trancheThreshold);
        int n = data.size();
        for ( int i = 0; i < n; i++ ) {
            if ( metric.getRunningMetric(i) >= metricThreshold ) {
                // we've found the largest group of variants with sensitivity >= our target truth sensitivity
                Tranche t = trancheOfVariants(data, i, trancheThreshold, model);
                logger.info(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, metric.getRunningMetric(i)));
                logger.info(String.format("  Tranche is %s", t));
                return t;
            }
        }

        return null;
    }

    public static Tranche trancheOfVariants(
            final List<VariantDatum> data,
            final int minI,
            final double ts,
            final VariantRecalibratorArgumentCollection.Mode model ) {
        int numKnown = 0, numNovel = 0, knownTi = 0, knownTv = 0, novelTi = 0, novelTv = 0;

        double minLod = data.get(minI).lod;
        for ( final VariantDatum datum : data ) {
            if ( datum.lod >= minLod ) {
                //if( ! datum.isKnown ) System.out.println(datum.pos);
                if ( datum.isKnown ) {
                    numKnown++;
                    if( datum.isSNP ) {
                        if ( datum.isTransition ) { knownTi++; } else { knownTv++; }
                    }
                } else {
                    numNovel++;
                    if( datum.isSNP ) {
                        if ( datum.isTransition ) { novelTi++; } else { novelTv++; }
                    }
                }
            }
        }

        double knownTiTv = knownTi / Math.max(1.0 * knownTv, 1.0);
        double novelTiTv = novelTi / Math.max(1.0 * novelTv, 1.0);

        int accessibleTruthSites = countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
        int nCallsAtTruth = countCallsAtTruth(data, minLod);

        return new Tranche(ts, minLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, nCallsAtTruth, model);
    }

    public static double fdrToTiTv(final double desiredFDR, final double targetTiTv) {
        return (1.0 - desiredFDR / 100.0) * (targetTiTv - 0.5) + 0.5;
    }

    public static int countCallsAtTruth(final List<VariantDatum> data, double minLOD ) {
        int n = 0;
        for ( VariantDatum d : data ) { n += (d.atTruthSite && d.lod >= minLOD ? 1 : 0); }
        return n;
    }
}
