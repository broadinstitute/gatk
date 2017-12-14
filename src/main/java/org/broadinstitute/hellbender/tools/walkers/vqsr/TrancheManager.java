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

        @Override
        public double getThreshold(final double tranche) {
            return 1.0 - tranche/100.0; // tranche of 1 => 99% sensitivity target
        }

        @Override
        public double getTarget() { return 1.0; }

        @Override
        public void calculateRunningMetric(final List<VariantDatum> data) {
            int nCalledAtTruth = 0;
            runningSensitivity = new double[data.size()];

            for ( int i = data.size() - 1; i >= 0; i-- ) {
                VariantDatum datum = data.get(i);
                nCalledAtTruth += datum.atTruthSite ? 1 : 0;
                runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
            }
        }

        @Override
        public double getRunningMetric(final int i) {
            return runningSensitivity[i];
        }

        @Override
        public int datumValue(final VariantDatum d) {
            return d.atTruthSite ? 1 : 0;
        }
    }

    public static List<TruthSensitivityTranche> findTranches(final List<VariantDatum> data,
                                                             final List<Double> tranches,
                                                             final SelectionMetric metric,
                                                             final VariantRecalibratorArgumentCollection.Mode model ) {
        return findTranches( data, tranches, metric, model, null );
    }

    public static List<TruthSensitivityTranche> findTranches(
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

        List<TruthSensitivityTranche> tranches = new ArrayList<>();
        for ( double trancheThreshold : trancheThresholds ) {
            TruthSensitivityTranche t = findTranche(data, metric, trancheThreshold, model);

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

    public static List<VQSLODTranche> findVQSLODTranches(
            final List<VariantDatum> data,
            final List<Double> trancheThresholds,
            final SelectionMetric metric,
            final VariantRecalibratorArgumentCollection.Mode model) {
        logger.info(String.format("Finding %d tranches for %d variants", trancheThresholds.size(), data.size()));

        Collections.sort( data, VariantDatum.VariantDatumLODComparator );
        metric.calculateRunningMetric(data);

        final List<VQSLODTranche> tranches = new ArrayList<>();
        for ( double trancheThreshold : trancheThresholds ) {
            VQSLODTranche t = findVQSLODTranche(data, metric, trancheThreshold, model);

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

    public static TruthSensitivityTranche findTranche(
            final List<VariantDatum> data,
            final SelectionMetric metric,
            final double trancheThreshold,
            final VariantRecalibratorArgumentCollection.Mode model ) {
        logger.info(String.format("  TruthSensitivityTranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, metric.getThreshold(trancheThreshold)));

        final double metricThreshold = metric.getThreshold(trancheThreshold);
        final int n = data.size();
        for ( int i = 0; i < n; i++ ) {
            if ( metric.getRunningMetric(i) >= metricThreshold ) {
                // we've found the largest group of variants with sensitivity >= our target truth sensitivity
                final TruthSensitivityTranche t = TruthSensitivityTranche.trancheOfVariants(data, i, trancheThreshold, model);
                logger.info(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, metric.getRunningMetric(i)));
                logger.info(String.format("  TruthSensitivityTranche is %s", t));
                return t;
            }
        }

        return null;
    }

    public static VQSLODTranche findVQSLODTranche(
            final List<VariantDatum> data,
            final SelectionMetric metric,
            final double trancheThreshold,
            final VariantRecalibratorArgumentCollection.Mode model ) {
        logger.info(String.format("  VQSLODTranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, metric.getThreshold(trancheThreshold)));

        final double metricThreshold = metric.getThreshold(trancheThreshold);
        final int n = data.size();
        for ( int i = 0; i < n; i++ ) {
            if ( data.get(i).lod >= trancheThreshold ) {
                // we've found the largest group of variants with LOD >= our target LOD
                final VQSLODTranche t = VQSLODTranche.trancheOfVariants(data, i, trancheThreshold, model);
                logger.info(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, metric.getRunningMetric(i)));
                logger.info(String.format("  VQSLODTranche is %s", t));
                return t;
            }
        }
        //If we don't have enough good data to meet the trancheThreshold, for VQSLODs we still want to return something
        logger.info(String.format("  Could not find tranche for %.3f: %.3f threshold; reporting empty tranche",
                trancheThreshold, metricThreshold));
        final VQSLODTranche t = VQSLODTranche.emptyTranche(data, n-1, trancheThreshold, model);
        return t;
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
