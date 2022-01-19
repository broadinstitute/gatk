package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

/**
 * TODO
 */
public class Tranche {

    static final String DEFAULT_TRANCHE_NAME = "anonymous";
    static final String COMMENT_STRING = "#";
    static final String VALUE_SEPARATOR = ",";
    static final int EXPECTED_COLUMN_COUNT = 9;

    protected final int accessibleTruthSites;
    protected final int callsAtTruthSites;
    final double minScore;  //minimum value of Score in this tranche
    final double novelTiTv;  //titv value of novel sites in this tranche
    final int numNovel;      //number of novel sites in this tranche
    final VariantTypeMode mode;
    final String name;       //Name of the tranche

    public Tranche(final String name,
                   final int numNovel,
                   final double minScore,
                   final VariantTypeMode mode,
                   final double novelTiTv,
                   final int accessibleTruthSites,
                   final int callsAtTruthSites) {
        if (numNovel < 0) {
            throw new GATKException("Invalid tranche - no. variants is < 0 : novel " + numNovel);
        }

        if (name == null) {
            throw new GATKException("BUG -- name cannot be null");
        }

        this.name = name;
        this.numNovel = numNovel;
        this.minScore = minScore;
        this.mode = mode;
        this.novelTiTv = novelTiTv;
        this.accessibleTruthSites = accessibleTruthSites;
        this.callsAtTruthSites = callsAtTruthSites;
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
     */
    static String tranchesString(final List<? extends Tranche> tranches) {
        try (final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
             final PrintStream stream = new PrintStream(bytes)) {
            if (tranches.size() > 1)
                tranches.sort(new TrancheComparator<>());

            Tranche prev = null;
            for (final Tranche t : tranches) {
                stream.print(t.getTrancheString(prev));
                prev = t;
            }

            return bytes.toString();
        }
        catch (final IOException e) {
            throw new GATKException("IOException while converting tranche to a string");
        }
    }

    public Double getTrancheIndex() {
        return getTruthSensitivity();
    }

    private <T extends Tranche> String getTrancheString(final T prev) {
            return String.format("%.2f,%d,%.4f,%.4f,VQSRTranche%s%.2fto%.2f,%s,%d,%d,%.4f%n",
                    getTrancheIndex(), numNovel, novelTiTv, minScore, mode.toString(),
                    (prev == null ? 0.0 : prev.getTrancheIndex()), getTrancheIndex(), mode.toString(), accessibleTruthSites, callsAtTruthSites, getTruthSensitivity());

    }

    static Tranche trancheOfVariants(final List<Double> scores,
                                     final List<Boolean> isTransition,
                                     final List<Boolean> isTruth,
                                     final int minI,
                                     final VariantTypeMode mode) {
        // TODO validate lengths
        int numNovel = 0, novelTi = 0, novelTv = 0;

        final double minScore = data.get(minI).score;
        for (final VariantDatum datum : data) {
            if (datum.score >= minScore) {
                numNovel++;
                if (datum.isTransition) {
                    novelTi++;
                } else {
                    novelTv++;
                }
            }
        }

        final double novelTiTv = novelTi / Math.max(1.0 * novelTv, 1.0);

        final int accessibleTruthSites = VariantDatum.countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
        final int nCallsAtTruth = VariantDatum.countCallsAtTruth(data, minScore);

        return new Tranche("unnamed", numNovel, minScore, mode, novelTiTv, accessibleTruthSites, nCallsAtTruth);
    }

    static double getRequiredDouble(final Map<String, String> bindings,
                                    final String key) {
        if (bindings.containsKey(key)) {
            try {
                return Double.parseDouble(bindings.get(key));
            } catch (final NumberFormatException e){
                throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
            }
        } else {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
    }

    static double getOptionalDouble(final Map<String, String> bindings, final String key, final double defaultValue) {
        try{
            return Double.parseDouble(bindings.getOrDefault(key, String.valueOf(defaultValue)));
        } catch (NumberFormatException e){
            throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
        }
    }

    static int getRequiredInteger(final Map<String, String> bindings, final String key) {
        if ( bindings.containsKey(key) ) {
            try{
                return Integer.parseInt(bindings.get(key));
            } catch (NumberFormatException e){
                throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
            }
        } else {
            throw new UserException.MalformedFile("Malformed tranches file.  Missing required key " + key);
        }
    }

    static int getOptionalInteger(final Map<String, String> bindings, final String key, final int defaultValue) {
        try{
            return Integer.parseInt(bindings.getOrDefault(key, String.valueOf(defaultValue)));
        } catch (NumberFormatException e){
            throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
        }
    }

    private double getTruthSensitivity() {
        return accessibleTruthSites > 0 ? callsAtTruthSites / (1.0*accessibleTruthSites) : 0.0;
    }
}
