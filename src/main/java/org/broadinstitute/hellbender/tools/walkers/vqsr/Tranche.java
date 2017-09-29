package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

/**
 * Created by gauthier on 7/13/17.
 */
public class Tranche {

    protected static final String DEFAULT_TRANCHE_NAME = "anonymous";
    protected static final String COMMENT_STRING = "#";
    protected static final String VALUE_SEPARATOR = ",";
    protected static final int EXPECTED_COLUMN_COUNT = 11;

    protected final int accessibleTruthSites;
    protected final int callsAtTruthSites;
    final double minVQSLod;  //minimum value of VQSLOD in this tranche
    final double knownTiTv;  //titv value of known sites in this tranche
    final double novelTiTv;  //titv value of novel sites in this tranche
    final int numKnown;      //number of known sites in this tranche
    final int numNovel;      //number of novel sites in this tranche
    final VariantRecalibratorArgumentCollection.Mode model;
    final String name;       //Name of the tranche

    public Tranche(final String name, final double knownTiTv, final int numNovel, final double minVQSLod, final VariantRecalibratorArgumentCollection.Mode model, final double novelTiTv, final int accessibleTruthSites, final int numKnown, final int callsAtTruthSites) {
        if ( numKnown < 0 || numNovel < 0) {
            throw new GATKException("Invalid tranche - no. variants is < 0 : known " + numKnown + " novel " + numNovel);
        }

        if ( name == null ) {
            throw new GATKException("BUG -- name cannot be null");
        }

        this.name = name;
        this.knownTiTv = knownTiTv;
        this.numNovel = numNovel;
        this.minVQSLod = minVQSLod;
        this.model = model;
        this.novelTiTv = novelTiTv;
        this.accessibleTruthSites = accessibleTruthSites;
        this.numKnown = numKnown;
        this.callsAtTruthSites = callsAtTruthSites;
    }

    public static class TrancheTruthSensitivityComparator implements Comparator<TruthSensitivityTranche> {
        @Override
        public int compare(final TruthSensitivityTranche tranche1, final TruthSensitivityTranche tranche2) {
            return Double.compare(tranche1.targetTruthSensitivity, tranche2.targetTruthSensitivity);
        }
    }

    public static class TrancheComparator<T extends Tranche> implements Comparator<T> {
        @Override
        public int compare(final T tranche1, final T tranche2) {
            //no matter what type of tranche we have, we want the output in order of increasing sensitivity, as measured by calls at truth sites
            return Double.compare(tranche1.callsAtTruthSites, tranche2.callsAtTruthSites);
        }
    }

    /**
     * Returns an appropriately formatted string representing the raw tranches file on disk.
     *
     * @param tranches
     * @return
     */
    public static String tranchesString( final List<? extends Tranche> tranches ) {
        try (final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
             final PrintStream stream = new PrintStream(bytes)) {
            if( tranches.size() > 1 )
                Collections.sort( tranches, new TrancheComparator<>());

            Tranche prev = null;
            for ( Tranche t : tranches ) {
                stream.printf(t.getTrancheString(prev));
                prev = t;
            }

            return bytes.toString();
        }
        catch (IOException e) {
            throw new GATKException("IOException while converting tranche to a string");
        }
    }

    public Double getTrancheIndex() {
        return getTruthSensitivity();
    }

    public <T extends Tranche> String getTrancheString(final T prev) {
            return String.format("%.2f,%d,%d,%.4f,%.4f,%.4f,VQSRTranche%s%.2fto%.2f,%s,%d,%d,%.4f%n",
                    getTrancheIndex(), numKnown, numNovel, knownTiTv, novelTiTv, minVQSLod, model.toString(),
                    (prev == null ? 0.0 : prev.getTrancheIndex()), getTrancheIndex(), model.toString(), accessibleTruthSites, callsAtTruthSites, getTruthSensitivity());

    }

    protected static Tranche trancheOfVariants(final List<VariantDatum> data, final int minI, final double ts, final VariantRecalibratorArgumentCollection.Mode model ) {
        int numKnown = 0, numNovel = 0, knownTi = 0, knownTv = 0, novelTi = 0, novelTv = 0;

        final double minLod = data.get(minI).lod;
        for (final VariantDatum datum : data) {
            if (datum.lod >= minLod) {
                if (datum.isKnown) {
                    numKnown++;
                    if (datum.isSNP) {
                        if (datum.isTransition) {
                            knownTi++;
                        } else {
                            knownTv++;
                        }
                    }
                } else {
                    numNovel++;
                    if (datum.isSNP) {
                        if (datum.isTransition) {
                            novelTi++;
                        } else {
                            novelTv++;
                        }
                    }
                }
            }
        }

        final double knownTiTv = knownTi / Math.max(1.0 * knownTv, 1.0);
        final double novelTiTv = novelTi / Math.max(1.0 * novelTv, 1.0);

        final int accessibleTruthSites = VariantDatum.countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
        final int nCallsAtTruth = VariantDatum.countCallsAtTruth(data, minLod);

        return new Tranche("unnamed", knownTiTv, numNovel, minLod, model, novelTiTv, accessibleTruthSites, numKnown, nCallsAtTruth);
    }

    protected static Tranche emptyTranche(final List<VariantDatum> data, final int minI, final double ts, final VariantRecalibratorArgumentCollection.Mode model ) {
        final double minLod = data.isEmpty() ? Double.NEGATIVE_INFINITY : data.get(minI).lod;
        final int accessibleTruthSites = VariantDatum.countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
        final int nCallsAtTruth = VariantDatum.countCallsAtTruth(data, minLod);

        final double knownTiTv = 0.0;
        final double novelTiTv = 0.0;
        final int numKnown = 0;
        final int numNovel = 0;

        return new Tranche("unnamed", knownTiTv, numNovel, minLod, model, novelTiTv, accessibleTruthSites, numKnown, nCallsAtTruth);
    }


    protected static double getRequiredDouble(final Map<String, String> bindings, final String key) {
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

    protected static double getOptionalDouble(final Map<String, String> bindings, String key, final double defaultValue) {
        try{
            return Double.valueOf(bindings.getOrDefault(key, String.valueOf(defaultValue)));
        } catch (NumberFormatException e){
            throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
        }
    }

    protected static int getRequiredInteger(final Map<String, String> bindings, final String key) {
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

    protected static int getOptionalInteger(final Map<String, String> bindings, final String key, final int defaultValue) {
        try{
            return Integer.valueOf(bindings.getOrDefault(key, String.valueOf(defaultValue)));
        } catch (NumberFormatException e){
            throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
        }
    }

    protected double getTruthSensitivity() {
        return accessibleTruthSites > 0 ? callsAtTruthSites / (1.0*accessibleTruthSites) : 0.0;
    }
}
