package org.broadinstitute.hellbender.tools.walkers.vqsr;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.File;
import java.io.IOException;
import java.util.*;

/*
 * Represents a truth sensitivity tranche in VQSR.
 * (Package-private because it's not usable outside.)
 */
final class Tranche {

    private static final String DEFAULT_TRANCHE_NAME = "anonymous";
    private static final String COMMENT_STRING = "#";
    private static final String VALUE_SEPARATOR = ",";
    private static final int EXPECTED_COLUMN_COUNT = 11;

    static final Comparator<Tranche> TRUTH_SENSITIVITY_ORDER = (tranche1, tranche2) -> Double.compare(tranche1.targetTruthSensitivity, tranche2.targetTruthSensitivity);

    private static final Logger logger = LogManager.getLogger(Tranche.class);

    //Note: visibility is set to package-local for testing
    final double targetTruthSensitivity;
    final double minVQSLod;  //minimum value of VQSLOD in this tranche
    final double knownTiTv;  //titv value of known sites in this tranche
    final double novelTiTv;  //titv value of novel sites in this tranche
    final int numKnown;      //number of known sites in this tranche
    final int numNovel;      //number of novel sites in this tranche
    final String name;       //Name of the tranche

    private final int accessibleTruthSites;
    private final int callsAtTruthSites;

    public Tranche(double targetTruthSensitivity, double minVQSLod, int numKnown, double knownTiTv, int numNovel, double novelTiTv, int accessibleTruthSites, int callsAtTruthSites, VariantRecalibratorArgumentCollection.Mode model, String name) {
        if ( targetTruthSensitivity < 0.0 || targetTruthSensitivity > 100.0) {
            throw new GATKException("Target FDR is unreasonable " + targetTruthSensitivity);
        }

        if ( numKnown < 0 || numNovel < 0) {
            throw new GATKException("Invalid tranche - no. variants is < 0 : known " + numKnown + " novel " + numNovel);
        }

        if ( name == null ) {
            throw new GATKException("BUG -- name cannot be null");
        }

        this.targetTruthSensitivity = targetTruthSensitivity;
        this.minVQSLod = minVQSLod;
        this.novelTiTv = novelTiTv;
        this.numNovel = numNovel;
        this.knownTiTv = knownTiTv;
        this.numKnown = numKnown;
        this.name = name;

        this.accessibleTruthSites = accessibleTruthSites;
        this.callsAtTruthSites = callsAtTruthSites;
    }

    @Override
    public String toString() {
        return String.format("Tranche targetTruthSensitivity=%.2f minVQSLod=%.4f known=(%d @ %.4f) novel=(%d @ %.4f) truthSites(%d accessible, %d called), name=%s]",
                targetTruthSensitivity, minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, name);
    }

    private static double getRequiredDouble(Map<String,String> bindings, String key ) {
        if ( bindings.containsKey(key) ) {
            try {
                return Double.valueOf(bindings.get(key));
            } catch (NumberFormatException e){
                throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
            }
        } else  {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
    }

    private static double getOptionalDouble(Map<String,String> bindings, String key, double defaultValue ) {
        try{
            return Double.valueOf(bindings.getOrDefault(key, String.valueOf(defaultValue)));
        } catch (NumberFormatException e){
            throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
        }
    }

    private static int getRequiredInteger(Map<String,String> bindings, String key) {
        if ( bindings.containsKey(key) ) {
            try{
                return Integer.valueOf(bindings.get(key));
            } catch (NumberFormatException e){
                throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
            }
        } else {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
    }

    private static int getOptionalInteger(Map<String,String> bindings, String key, int defaultValue) {
        try{
            return Integer.valueOf(bindings.getOrDefault(key, String.valueOf(defaultValue)));
        } catch (NumberFormatException e){
            throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
        }
    }

    /**
     * Returns a list of tranches, sorted from most to least specific, read in from file f.
     * @throws IOException if there are problems reading the file.
     */
    public static List<Tranche> readTranches(File f) throws IOException{
        String[] header = null;
        List<Tranche> tranches = new ArrayList<>();

        try (XReadLines xrl = new XReadLines(f) ) {
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

                    Map<String, String> bindings = new HashMap<>();
                    for (int i = 0; i < vals.length; i++) {
                        bindings.put(header[i], vals[i]);
                    }
                    tranches.add(new Tranche(
                            getRequiredDouble(bindings, "targetTruthSensitivity"),
                            getRequiredDouble(bindings, "minVQSLod"),
                            getOptionalInteger(bindings, "numKnown", -1),
                            getOptionalDouble(bindings, "knownTiTv", -1.0),
                            getRequiredInteger(bindings, "numNovel"),
                            getRequiredDouble(bindings, "novelTiTv"),
                            getOptionalInteger(bindings, "accessibleTruthSites", -1),
                            getOptionalInteger(bindings, "callsAtTruthSites", -1),
                            VariantRecalibratorArgumentCollection.Mode.valueOf(bindings.get("model")),
                            bindings.get("filterName")));
                }
            }

            tranches.sort(TRUTH_SENSITIVITY_ORDER);
            return tranches;
        }
    }

    @VisibleForTesting
    static final class TruthSensitivityMetric {
        private final String name;
        private double[] runningSensitivity;
        private final int nTrueSites;

        public TruthSensitivityMetric(int nTrueSites) {
            this.name = "TruthSensitivity";
            this.nTrueSites = nTrueSites;
        }

        public String getName(){
            return name;
        }

        public double getThreshold(double tranche) {
            return 1.0 - tranche/100.0; // tranche of 1 => 99% sensitivity target
        }

        public void calculateRunningMetric(List<VariantDatum> data) {
            int nCalledAtTruth = 0;
            runningSensitivity = new double[data.size()];

            for ( int i = data.size() - 1; i >= 0; i-- ) {
                VariantDatum datum = data.get(i);
                nCalledAtTruth += datum.atTruthSite ? 1 : 0;
                runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
            }
        }

        public double getRunningMetric(int i) {
            return runningSensitivity[i];
        }
    }

    public static List<Tranche> findTranches( final List<VariantDatum> data, final double[] trancheThresholds, final TruthSensitivityMetric metric, final VariantRecalibratorArgumentCollection.Mode model) {
        logger.info(String.format("Finding %d tranches for %d variants", trancheThresholds.length, data.size()));

        Collections.sort(data, VariantDatum.VariantDatumLODComparator);
        metric.calculateRunningMetric(data);

        List<Tranche> tranches = new ArrayList<>();
        for ( double trancheThreshold : trancheThresholds ) {
            Tranche t = findTranche(data, metric, trancheThreshold, model);

            if ( t == null ) {
                if ( tranches.size() == 0 ) {
                    throw new UserException(String.format("Couldn't find any tranche containing variants with a %s > %.2f. Are you sure the truth files contain unfiltered variants which overlap the input data?", metric.getName(), metric.getThreshold(trancheThreshold)));
                }
                break;
            }

            tranches.add(t);
        }

        tranches.sort(Tranche.TRUTH_SENSITIVITY_ORDER);
        return tranches;
    }

    private static Tranche findTranche( final List<VariantDatum> data, final TruthSensitivityMetric metric, final double trancheThreshold, final VariantRecalibratorArgumentCollection.Mode model ) {
        logger.debug(String.format("  Tranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, metric.getThreshold(trancheThreshold)));

        double metricThreshold = metric.getThreshold(trancheThreshold);
        int n = data.size();
        for ( int i = 0; i < n; i++ ) {
            if ( metric.getRunningMetric(i) >= metricThreshold ) {
                // we've found the largest group of variants with sensitivity >= our target truth sensitivity
                Tranche t = trancheOfVariants(data, i, trancheThreshold, model);
                logger.debug(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, metric.getRunningMetric(i)));
                logger.debug(String.format("  Tranche is %s", t));
                return t;
            }
        }

        return null;
    }

    private static Tranche trancheOfVariants( final List<VariantDatum> data, int minI, double ts, final VariantRecalibratorArgumentCollection.Mode model ) {
        int numKnown = 0, numNovel = 0, knownTi = 0, knownTv = 0, novelTi = 0, novelTv = 0;

        double minLod = data.get(minI).lod;
        for ( final VariantDatum datum : data ) {
            if ( datum.lod >= minLod ) {
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

        int accessibleTruthSites = VariantDatum.countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
        int nCallsAtTruth = VariantDatum.countCallsAtTruth(data, minLod);

        return new Tranche(ts, minLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, nCallsAtTruth, model, DEFAULT_TRANCHE_NAME);
    }

}
