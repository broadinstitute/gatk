package org.broadinstitute.hellbender.tools.walkers.vqsr;

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

    private final static Logger logger = LogManager.getLogger(Tranche.class);

    //Note: visibility is set to package-local for testing
    final double targetTruthSensitivity;
    final double minVQSLod;  //minimum value of VQSLOD in this tranche
    final double knownTiTv;  //titv value of known sites in this tranche
    final double novelTiTv;  //titv value of novel sites in this tranche
    final int numKnown;      //number of known sites in this tranche
    final int numNovel;      //number of novel sites in this tranche
    final String name;       //Name of the tranche
    private final VariantRecalibratorArgumentCollection.Mode model;    //this is a SNP VQSR tranche or Indel tranche?

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
        this.model = model;
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
     * @throws java.io.IOException if there are problems reading the file.
     */
    public static List<Tranche> readTranches(File f) throws IOException{
        String[] header = null;
        final List<Tranche> tranches = new ArrayList<>();

        try (XReadLines xrl = new XReadLines(f) ) {
            for (final String line : xrl) {
                if (line.startsWith(COMMENT_STRING)) {
                    continue;
                }

                final String[] vals = line.split(VALUE_SEPARATOR);
                if (header == null) {  //reading the header
                    header = vals;
                    if (header.length != EXPECTED_COLUMN_COUNT) {
                        throw new UserException.MalformedFile(f, "Expected " + EXPECTED_COLUMN_COUNT + " elements in header line " + line);
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

    /**
     * Creates and returns a list of variant tranches for given levels
     * of thresholds in truth sensitivity.
     */
    public static List<Tranche> findTranches( final List<VariantDatum> data, final double[] trancheThresholds, final int nCallsAtTruth, final VariantRecalibratorArgumentCollection.Mode model) {
        logger.info(String.format("Finding %d tranches for %d variants", trancheThresholds.length, data.size()));

        Collections.sort(data, VariantDatum.VariantDatumLODComparator);
        final double[] runningSensitivity = calculateRunningSensitivity(data, nCallsAtTruth);

        final List<Tranche> tranches = new ArrayList<>(trancheThresholds.length);
        for ( double trancheThreshold : trancheThresholds ) {
            final Tranche t = findTranche(data, runningSensitivity, trancheThreshold, model);

            if ( t == null ) {
                if ( tranches.isEmpty() ) {
                    throw new UserException(String.format("Couldn't find any tranche containing variants with TruthSensitivity > %.2f. " +
                            "Are you sure the truth files contain unfiltered variants which overlap the input data?", getSensitivityThreshold(trancheThreshold)));
                }
                break;
            }

            tranches.add(t);
        }

        tranches.sort(Tranche.TRUTH_SENSITIVITY_ORDER);
        return tranches;
    }

    public static double getSensitivityThreshold(double tranche) {
        return 1.0 - tranche/100.0; // tranche of 1 => 99% sensitivity target
    }

    /**
     * Given a list of data points sorted by quality (we use LOD), returns the array
     * of sensitivity scores of variants with LOD scored above the given.
     * That is, the resulting array has the same length as the input list and at each position i
     * it contains the sensitivity of variants with indices at least equal to i.
     */
    private static double[] calculateRunningSensitivity(final List<VariantDatum> data, final int nTrueSites) {
        final double[] runningSensitivity = new double[data.size()];
        int nCalledAtTruth = 0;

        for ( int i = data.size() - 1; i >= 0; i-- ) {
            VariantDatum datum = data.get(i);
            nCalledAtTruth += datum.atTruthSite ? 1 : 0;
            runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
        }
        return runningSensitivity;
    }

    /**
     * Creates a tranche of variants that have the truth sensitivity better than the threshold.
     * Returns null if no such tranche exists.
     */
    private static Tranche findTranche( final List<VariantDatum> data, final double[] runningSensitivity, final double trancheThreshold, final VariantRecalibratorArgumentCollection.Mode model ) {
        final double metricThreshold = getSensitivityThreshold(trancheThreshold);
        logger.debug(String.format("  Tranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, metricThreshold));

        for ( int i = 0; i < data.size(); i++ ) {
            if ( runningSensitivity[i] >= metricThreshold ) {
                // we've found the largest group of variants with sensitivity >= our target truth sensitivity
                Tranche t = trancheOfVariants(data, data.get(i).lod, trancheThreshold, model);
                logger.debug(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, runningSensitivity[i]));
                logger.debug(String.format("  Tranche is %s", t));
                return t;
            }
        }

        return null;
    }

    /**
     * Creates the tranche of variants with with minimal LOD given by {@param minLod}.
     */
    private static Tranche trancheOfVariants( final List<VariantDatum> data, final double minLod, final double ts, final VariantRecalibratorArgumentCollection.Mode model ) {
        int numKnown = 0, numNovel = 0, knownTi = 0, knownTv = 0, novelTi = 0, novelTv = 0;

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

        final double knownTiTv = knownTi / Math.max(1.0 * knownTv, 1.0);
        final double novelTiTv = novelTi / Math.max(1.0 * novelTv, 1.0);

        final int accessibleTruthSites = VariantDatum.countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
        final int nCallsAtTruth = VariantDatum.countCallsAtTruth(data, minLod);

        return new Tranche(ts, minLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, nCallsAtTruth, model, DEFAULT_TRANCHE_NAME);
    }

}
