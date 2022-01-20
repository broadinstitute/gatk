package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/*
 * TODO this whole class needs heavy refactoring. it's cleaned up significantly from the multiple tranche classes in VQSR,
 *  but still has a long way to go. we should decide how strongly to couple to VariantDatum/VariantDataManger and refactor
 *  those classes at the same time
 */
final class TruthSensitivityTranche {
    private static final int CURRENT_VERSION = 6;

    private static final Comparator<TruthSensitivityTranche> TRUTH_SENSITIVITY_ORDER = Comparator.comparingDouble(tranche -> tranche.targetTruthSensitivity);

    private static final Logger logger = LogManager.getLogger(TruthSensitivityTranche.class);

    private static final String DEFAULT_TRANCHE_NAME = "anonymous";
    private static final String COMMENT_STRING = "#";
    private static final String VALUE_SEPARATOR = ",";
    private static final int EXPECTED_COLUMN_COUNT = 9;

    private final String name;
    private final int numSites;
    private final double minScore;
    private final VariantTypeMode mode;
    private final double tiTv;
    private final int accessibleTruthSites;
    private final int callsAtTruthSites;
    private final double targetTruthSensitivity;

    private TruthSensitivityTranche(final double targetTruthSensitivity,
                                    final int numSites,
                                    final double tiTv,
                                    final double minScore,
                                    final String name,
                                    final VariantTypeMode mode,
                                    final int accessibleTruthSites,
                                    final int callsAtTruthSites) {
        // TODO more validation
        Utils.nonNull(name);
        ParamUtils.isPositiveOrZero(numSites, "Number of variants cannot be negative.");
        ParamUtils.inRange(targetTruthSensitivity, 0., 100.,"Target truth sensitivity must be in [0, 100].");

        this.targetTruthSensitivity = targetTruthSensitivity;
        this.numSites = numSites;
        this.tiTv = tiTv;
        this.minScore = minScore;
        this.name = name;
        this.mode = mode;
        this.accessibleTruthSites = accessibleTruthSites;
        this.callsAtTruthSites = callsAtTruthSites;
    }

    public static String printHeader() {
        try (final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
             final PrintStream stream = new PrintStream(bytes)) {

            stream.println("# Variant quality score tranches file");
            stream.println("# Version number " + CURRENT_VERSION);
            stream.println("targetTruthSensitivity,numSites,tiTv,minScore,filterName,mode,accessibleTruthSites,callsAtTruthSites,truthSensitivity");

            return bytes.toString();
        }
        catch (IOException e) {
            throw new GATKException("IOException while converting tranche to a string");
        }
    }

    @Override
    public String toString() {
        return String.format("TruthSensitivityTranche targetTruthSensitivity=%.2f minScore=%.4f numSites=(%d @ %.4f) truthSites(%d accessible, %d called), name=%s]",
                targetTruthSensitivity, minScore, numSites, tiTv, accessibleTruthSites, callsAtTruthSites, name);
    }

    /**
     * Returns a list of tranches, sorted from most to least specific, read in from file f.
     * @throws IOException if there are problems reading the file.
     */
    static List<TruthSensitivityTranche> readTranches(final GATKPath f) throws IOException{
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
                        throw new UserException.MalformedFile(f, String.format("Expected %d elements in header line %s", EXPECTED_COLUMN_COUNT, line));
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
                            getRequiredInteger(bindings, "numSites"),
                            getRequiredDouble(bindings, "tiTv"),
                            getRequiredDouble(bindings, "minScore"),
                            bindings.get("filterName"),
                            VariantTypeMode.valueOf(bindings.get("mode")),
                            getOptionalInteger(bindings, "accessibleTruthSites", -1),
                            getOptionalInteger(bindings, "callsAtTruthSites", -1)
                    ));
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

        TruthSensitivityMetric(final int nTrueSites) {
            this.name = "TruthSensitivity";
            this.nTrueSites = nTrueSites;
        }

        public String getName(){
            return name;
        }

        double getRunningMetric(final int i) {
            return runningSensitivity[i];
        }

        private static double getThreshold(final double tranche) {
            return 1. - tranche / 100.; // tranche of 1 => 99% sensitivity target
        }

        private void calculateRunningMetric(final List<Boolean> isTruth) {
            int nCalledAtTruth = 0;
            runningSensitivity = new double[isTruth.size()];

            for (int i = isTruth.size() - 1; i >= 0; i--) {
                nCalledAtTruth += isTruth.get(i) ? 1 : 0;
                runningSensitivity[i] = 1 - nCalledAtTruth / (1. * nTrueSites);
            }
        }
    }

    // TODO clean all this up once VariantDataManager is refactored
    static List<TruthSensitivityTranche> findTranches(final List<Double> scores,
                                                      final List<Boolean> isBiallelicSNP,
                                                      final List<Boolean> isTransition,
                                                      final List<Boolean> isTruth,
                                                      final List<Double> trancheThresholds,
                                                      final TruthSensitivityMetric metric,
                                                      final VariantTypeMode mode) {
        // TODO validate lengths
        logger.info(String.format("Finding %d tranches for %d variants", trancheThresholds.size(), scores.size()));

        final List<Integer> indicesSortedByScore = IntStream.range(0, scores.size()).boxed()
                .sorted(Comparator.comparingDouble(scores::get))
                .collect(Collectors.toList());
        final List<Double> sortedScores = indicesSortedByScore.stream().map(scores::get).collect(Collectors.toList());
        final List<Boolean> sortedIsBiallelicSNP = indicesSortedByScore.stream().map(isBiallelicSNP::get).collect(Collectors.toList());
        final List<Boolean> sortedIsTransition = indicesSortedByScore.stream().map(isTransition::get).collect(Collectors.toList());
        final List<Boolean> sortedIsTruth = indicesSortedByScore.stream().map(isTruth::get).collect(Collectors.toList());
        metric.calculateRunningMetric(sortedIsTruth);

        List<TruthSensitivityTranche> tranches = new ArrayList<>(trancheThresholds.size());
        for (double trancheThreshold : trancheThresholds) {
            TruthSensitivityTranche t = findTranche(sortedScores, sortedIsBiallelicSNP, sortedIsTransition, sortedIsTruth, metric, trancheThreshold, mode);

            if ( t == null ) {
                if (tranches.isEmpty()) {
                    throw new UserException(String.format("Couldn't find any tranche containing variants with a %s > %.2f. Are you sure the truth files contain unfiltered variants which overlap the input data?", metric.getName(), TruthSensitivityMetric.getThreshold(trancheThreshold)));
                }
                break;
            }

            tranches.add(t);
        }

        tranches.sort(TRUTH_SENSITIVITY_ORDER);
        return tranches;
    }

    private static TruthSensitivityTranche findTranche(final List<Double> sortedScores,
                                                       final List<Boolean> sortedIsBiallelicSNP,
                                                       final List<Boolean> sortedIsTransition,
                                                       final List<Boolean> sortedIsTruth,
                                                       final TruthSensitivityMetric metric,
                                                       final double trancheThreshold,
                                                       final VariantTypeMode mode) {
        // TODO validate lengths
        logger.debug(String.format("  TruthSensitivityTranche threshold %.2f => selection metric threshold %.3f", trancheThreshold, TruthSensitivityMetric.getThreshold(trancheThreshold)));

        double metricThreshold = TruthSensitivityMetric.getThreshold(trancheThreshold);
        int n = sortedScores.size();
        for (int i = 0; i < n; i++) {
            if (metric.getRunningMetric(i) >= metricThreshold) {
                // we've found the largest group of variants with sensitivity >= our target truth sensitivity
                TruthSensitivityTranche t = trancheOfVariants(sortedScores, sortedIsBiallelicSNP, sortedIsTransition, sortedIsTruth, i, trancheThreshold, mode);
                logger.debug(String.format("  Found tranche for %.3f: %.3f threshold starting with variant %d; running score is %.3f ",
                        trancheThreshold, metricThreshold, i, metric.getRunningMetric(i)));
                logger.debug(String.format("  TruthSensitivityTranche is %s", t));
                return t;
            }
        }

        return null;
    }

    private static TruthSensitivityTranche trancheOfVariants(final List<Double> sortedScores,
                                                             final List<Boolean> sortedIsBiallelicSNP,
                                                             final List<Boolean> sortedIsTransition,
                                                             final List<Boolean> sortedIsTruth,
                                                             final int minI,
                                                             final double targetTruthSensitivity,
                                                             final VariantTypeMode mode) {
        int numSites = 0;
        int ti = 0;
        int tv = 0;

        final double minScore = sortedScores.get(minI);
        for (int i = 0; i < sortedScores.size(); i++) {
            if (sortedScores.get(i) >= minScore) {
                numSites++;
                if (sortedIsBiallelicSNP.get(i)) {
                    if (sortedIsTransition.get(i)) {
                        ti++;
                    } else {
                        tv++;
                    }
                }
            }
        }

        final double tiTv = ti / Math.max(tv, 1.);

        final int accessibleTruthSites = countCallsAtTruth(sortedScores, sortedIsTruth, Double.NEGATIVE_INFINITY);
        final int callsAtTruthSites = countCallsAtTruth(sortedScores, sortedIsTruth, minScore);

        return new TruthSensitivityTranche(targetTruthSensitivity, numSites, tiTv, minScore, DEFAULT_TRANCHE_NAME, mode, accessibleTruthSites, callsAtTruthSites);
    }

    public static class TrancheComparator implements Comparator<TruthSensitivityTranche> {
        @Override
        public int compare(final TruthSensitivityTranche tranche1,
                           final TruthSensitivityTranche tranche2) {
            //no matter what type of tranche we have, we want the output in order of increasing sensitivity, as measured by calls at truth sites
            return Double.compare(tranche1.callsAtTruthSites, tranche2.callsAtTruthSites);
        }
    }

    /**
     * Returns an appropriately formatted string representing the raw tranches file on disk.
     */
    static String tranchesString(final List<TruthSensitivityTranche> tranches) {
        try (final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
             final PrintStream stream = new PrintStream(bytes)) {
            if (tranches.size() > 1)
                tranches.sort(new TrancheComparator());

            TruthSensitivityTranche prev = null;
            for (final TruthSensitivityTranche t : tranches) {
                stream.print(t.getTrancheString(prev));
                prev = t;
            }

            return bytes.toString();
        }
        catch (final IOException e) {
            throw new GATKException("IOException while converting tranche to a string");
        }
    }

    private String getTrancheString(final TruthSensitivityTranche prev) {
        return String.format("%.2f,%d,%.4f,%.4f,VQSRTranche%s%.2fto%.2f,%s,%d,%d,%.4f%n",
                targetTruthSensitivity, numSites, tiTv, minScore, mode.toString(),
                (prev == null ? 0.0 : prev.targetTruthSensitivity), targetTruthSensitivity, mode.toString(), accessibleTruthSites, callsAtTruthSites, getTruthSensitivity());

    }

    private static double getRequiredDouble(final Map<String, String> bindings,
                                            final String key) {
        if (bindings.containsKey(key)) {
            try {
                return Double.parseDouble(bindings.get(key));
            } catch (final NumberFormatException e){
                throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
            }
        } else {
            throw new UserException.MalformedFile("Malformed tranches file. Missing required key " + key);
        }
    }

    private static int getRequiredInteger(final Map<String, String> bindings,
                                          final String key) {
        if (bindings.containsKey(key)) {
            try{
                return Integer.parseInt(bindings.get(key));
            } catch (final NumberFormatException e){
                throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
            }
        } else {
            throw new UserException.MalformedFile("Malformed tranches file. Missing required key " + key);
        }
    }

    private static int getOptionalInteger(final Map<String, String> bindings, final String key,
                                          final int defaultValue) {
        try {
            return Integer.parseInt(bindings.getOrDefault(key, String.valueOf(defaultValue)));
        } catch (final NumberFormatException e){
            throw new UserException.MalformedFile("Malformed tranches file. Invalid value for key " + key);
        }
    }

    private double getTruthSensitivity() {
        return accessibleTruthSites > 0 ? callsAtTruthSites / ((double) accessibleTruthSites) : 0.;
    }

    static int countCallsAtTruth(final List<Double> scores,
                                 final List<Boolean> isTruth,
                                 final double minScore) {
        return (int) IntStream.range(0, scores.size()).filter(i -> isTruth.get(i) && scores.get(i) >= minScore).count();
    }

    public String getName() {
        return name;
    }

    public int getNumSites() {
        return numSites;
    }

    public double getMinScore() {
        return minScore;
    }

    public VariantTypeMode getMode() {
        return mode;
    }

    public double getTiTv() {
        return tiTv;
    }

    public int getAccessibleTruthSites() {
        return accessibleTruthSites;
    }

    public int getCallsAtTruthSites() {
        return callsAtTruthSites;
    }

    public double getTargetTruthSensitivity() {
        return targetTruthSensitivity;
    }
}
