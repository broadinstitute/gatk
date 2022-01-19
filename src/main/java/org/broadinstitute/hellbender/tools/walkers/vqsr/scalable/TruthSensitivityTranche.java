package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/*
 * TODO
 */
final class TruthSensitivityTranche extends Tranche {
    private static final int CURRENT_VERSION = 5;

    private static final Comparator<TruthSensitivityTranche> TRUTH_SENSITIVITY_ORDER = (tranche1, tranche2) -> Double.compare(tranche1.targetTruthSensitivity, tranche2.targetTruthSensitivity);

    private static final Logger logger = LogManager.getLogger(TruthSensitivityTranche.class);

    //Note: visibility is set to package-local for testing
    final double targetTruthSensitivity;

    public TruthSensitivityTranche(
            final double targetTruthSensitivity,
            final double minScore,
            final int numNovel,
            final double novelTiTv,
            final int accessibleTruthSites,
            final int callsAtTruthSites,
            final VariantTypeMode mode) {
        this(targetTruthSensitivity, minScore, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, mode, "anonymous");
    }

    public TruthSensitivityTranche(
            final double targetTruthSensitivity,
            final double minScore,
            final int numNovel,
            final double novelTiTv,
            final int accessibleTruthSites,
            final int callsAtTruthSites,
            final VariantTypeMode mode,
            final String name) {
        super(name, numNovel, minScore, mode, novelTiTv, accessibleTruthSites, callsAtTruthSites);
        if (targetTruthSensitivity < 0.0 || targetTruthSensitivity > 100.0) {
            throw new GATKException("Target FDR is unreasonable " + targetTruthSensitivity);
        }

        if (numNovel < 0) {
            throw new GATKException("Invalid tranche - no. variants is < 0 : novel " + numNovel);
        }

        if (name == null) {
            throw new GATKException("BUG -- name cannot be null");
        }

        this.targetTruthSensitivity = targetTruthSensitivity;

    }

    public Double getTrancheIndex() {
        return targetTruthSensitivity;
    }

    public static String printHeader() {
        try (final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
             final PrintStream stream = new PrintStream(bytes)) {

            stream.println("# Variant quality score tranches file");
            stream.println("# Version number " + CURRENT_VERSION);
            stream.println("targetTruthSensitivity,numNovel,novelTiTv,minScore,filterName,mode,accessibleTruthSites,callsAtTruthSites,truthSensitivity");

            return bytes.toString();
        }
        catch (IOException e) {
            throw new GATKException("IOException while converting tranche to a string");
        }
    }

    @Override
    public String toString() {
        return String.format("TruthSensitivityTranche targetTruthSensitivity=%.2f minScore=%.4f novel=(%d @ %.4f) truthSites(%d accessible, %d called), name=%s]",
                targetTruthSensitivity, minScore, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, name);
    }

    /**
     * Returns a list of tranches, sorted from most to least specific, read in from file f.
     * @throws IOException if there are problems reading the file.
     */
    public static List<TruthSensitivityTranche> readTranches(final GATKPath f) throws IOException{
        String[] header = null;
        final List<TruthSensitivityTranche> tranches = new ArrayList<>();

        try (XReadLines xrl = new XReadLines(f.toPath())) {
            for (final String line : xrl) {
                if (line.startsWith(COMMENT_STRING)) {
                    continue;
                }

                final String[] vals = line.split(VALUE_SEPARATOR);
                if (header == null) {  //reading the header
                    header = vals;
                    if (header.length != EXPECTED_COLUMN_COUNT) {
                        throw new UserException.MalformedFile(f, "Expected 11 elements in header line " + line);
                    }
                } else {
                    if (header.length != vals.length) {
                        throw new UserException.MalformedFile(f, "Line had too few/many fields.  Header = " + header.length + " vals " + vals.length + ". The line was: " + line);
                    }

                    Map<String, String> bindings = new LinkedHashMap<>();
                    for (int i = 0; i < vals.length; i++) {
                        bindings.put(header[i], vals[i]);
                    }
                    tranches.add(new TruthSensitivityTranche(
                            getRequiredDouble(bindings, "targetTruthSensitivity"),
                            getRequiredDouble(bindings, "minScore"),
                            getRequiredInteger(bindings, "numNovel"),
                            getRequiredDouble(bindings, "novelTiTv"),
                            getOptionalInteger(bindings, "accessibleTruthSites", -1),
                            getOptionalInteger(bindings, "callsAtTruthSites", -1),
                            VariantTypeMode.valueOf(bindings.get("mode")),
                            bindings.get("filterName")));
                }
            }
        }

        tranches.sort(TRUTH_SENSITIVITY_ORDER);
        return tranches;
    }

    @VisibleForTesting
    static final class TruthSensitivityMetric {
        private final String name;
        private double[] runningSensitivity;
        private final int nTrueSites;

        public TruthSensitivityMetric(final int nTrueSites) {
            this.name = "TruthSensitivity";
            this.nTrueSites = nTrueSites;
        }

        public String getName(){
            return name;
        }

        public double getThreshold(final double tranche) {
            return 1.0 - tranche/100.0; // tranche of 1 => 99% sensitivity target
        }

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
    }

    public static List<TruthSensitivityTranche> findTranches(final List<VariantDatum> data, final double[] trancheThresholds, final TruthSensitivityMetric metric, final VariantTypeMode mode) {
        logger.info(String.format("Finding %d tranches for %d variants", trancheThresholds.length, data.size()));

        data.sort(VariantDatum.VariantDatumLODComparator);
        metric.calculateRunningMetric(data);

        List<TruthSensitivityTranche> tranches = new ArrayList<>();
        for ( double trancheThreshold : trancheThresholds ) {
            TruthSensitivityTranche t = findTranche(data, metric, trancheThreshold, mode);

            if ( t == null ) {
                if (tranches.isEmpty()) {
                    throw new UserException(String.format("Couldn't find any tranche containing variants with a %s > %.2f. Are you sure the truth files contain unfiltered variants which overlap the input data?", metric.getName(), metric.getThreshold(trancheThreshold)));
                }
                break;
            }

            tranches.add(t);
        }

        tranches.sort(TruthSensitivityTranche.TRUTH_SENSITIVITY_ORDER);
        return tranches;
    }

    private static TruthSensitivityTranche findTranche(final List<Double> scores,
                                                       final List<Boolean> isTransition,
                                                       final List<Boolean> isTruth,
                                                       final TruthSensitivityMetric metric,
                                                       final double trancheThreshold,
                                                       final VariantTypeMode mode ) {
        logger.debug(String.format("  TruthSensitivityTranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, metric.getThreshold(trancheThreshold)));

        double metricThreshold = metric.getThreshold(trancheThreshold);
        int n = data.size();
        for ( int i = 0; i < n; i++ ) {
            if ( metric.getRunningMetric(i) >= metricThreshold ) {
                // we've found the largest group of variants with sensitivity >= our target truth sensitivity
                TruthSensitivityTranche t = trancheOfVariants(scores, isTransition, isTruth, i, trancheThreshold, mode);
                logger.debug(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, metric.getRunningMetric(i)));
                logger.debug(String.format("  TruthSensitivityTranche is %s", t));
                return t;
            }
        }

        return null;
    }

    protected static TruthSensitivityTranche trancheOfVariants(final List<Double> scores,
                                                               final List<Boolean> isTransition,
                                                               final List<Boolean> isTruth,
                                                               final int minI,
                                                               final double ts,
                                                               final VariantTypeMode mode ) {
        final Tranche basicTranche = Tranche.trancheOfVariants(scores, isTransition, isTruth, minI, mode);
        return new TruthSensitivityTranche(ts, basicTranche.minScore, basicTranche.numNovel, basicTranche.novelTiTv, basicTranche.accessibleTruthSites, basicTranche.callsAtTruthSites, mode, DEFAULT_TRANCHE_NAME);
    }

}
