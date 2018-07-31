package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import static java.lang.Math.abs;
import static java.util.Collections.sort;

/*
 * Represents a VQSLOD tranche in VQSR for use in scattered VariantRecalibrator runs.
 * (Package-private because it's not usable outside.)
 */
public class VQSLODTranche extends Tranche {
    private static final int CURRENT_VERSION = 6;

    public Double getTrancheIndex() {
        return minVQSLod;
    }

    public VQSLODTranche(
            final double minVQSLod,
            final int numKnown,
            final double knownTiTv,
            final int numNovel,
            final double novelTiTv,
            final int accessibleTruthSites,
            final int callsAtTruthSites,
            final VariantRecalibratorArgumentCollection.Mode model,
            final String name) {
        super(name, knownTiTv, numNovel, minVQSLod, model, novelTiTv, accessibleTruthSites, numKnown, callsAtTruthSites);
    }

    public int compareTo(VQSLODTranche t) {
        return (int)Math.round(t.minVQSLod - this.minVQSLod);
    }

    @Override
    public String toString() {
        return String.format("Tranche minVQSLod=%.4f known=(%d @ %.4f) novel=(%d @ %.4f) truthSites(%d accessible, %d called), name=%s]",
                minVQSLod, numKnown, knownTiTv, numNovel, novelTiTv, accessibleTruthSites, callsAtTruthSites, name);
    }

    public static String printHeader() {
        try (final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
             final PrintStream stream = new PrintStream(bytes)) {

            stream.println("# Variant quality score tranches file");
            stream.println("# Version number " + CURRENT_VERSION);
            stream.println("requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity");

            return bytes.toString();
        }
        catch (IOException e) {
            throw new GATKException("IOException while converting tranche to a string");
        }
    }

    protected static VQSLODTranche trancheOfVariants(final List<VariantDatum> data, final int minI, final double trancheThreshold, final VariantRecalibratorArgumentCollection.Mode model ) {
        final Tranche basicTranche = Tranche.trancheOfVariants(data, minI, trancheThreshold, model);

        //First column should be the requested threshold, not the value in the data closest to the threshold
        return new VQSLODTranche(trancheThreshold, basicTranche.numKnown, basicTranche.knownTiTv, basicTranche.numNovel, basicTranche.novelTiTv, basicTranche.accessibleTruthSites, basicTranche.callsAtTruthSites, model, DEFAULT_TRANCHE_NAME);
    }

    protected static VQSLODTranche emptyTranche(final List<VariantDatum> data, final int minI, final double trancheThreshold, final VariantRecalibratorArgumentCollection.Mode model ) {
        final Tranche basicTranche = Tranche.emptyTranche(data, minI, trancheThreshold, model);

        //First column should be the requested threshold, not the value in the data closest to the threshold
        return new VQSLODTranche(trancheThreshold, basicTranche.numKnown, basicTranche.knownTiTv, basicTranche.numNovel, basicTranche.novelTiTv, basicTranche.accessibleTruthSites, basicTranche.callsAtTruthSites, model, DEFAULT_TRANCHE_NAME);
    }

    /**
     * Returns a list of tranches, sorted from most to least specific, read in from file f.
     * @throws IOException if there are problems reading the file.
     */
    public static List<VQSLODTranche> readTranches(final File f) throws IOException{
        String[] header = null;
        List<VQSLODTranche> tranches = new ArrayList<>();

        try (XReadLines xrl = new XReadLines(f) ) {
            for (final String line : xrl) {
                if (line.startsWith(COMMENT_STRING)) {
                    if ( !line.contains("Version"))
                        continue;
                    else {
                        String[] words = line.split("\\s+"); //split on whitespace
                        if (Integer.parseInt(words[3]) != CURRENT_VERSION)
                            throw new UserException.BadInput("The file " + " contains version " + words[3] + " tranches, but VQSLOD tranche parsing requires version " + CURRENT_VERSION);
                        continue;
                    }
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

                    Map<String, String> bindings = new LinkedHashMap<>();
                    for (int i = 0; i < vals.length; i++) {
                        bindings.put(header[i], vals[i]);
                    }
                    tranches.add(new VQSLODTranche(
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
        }

        tranches.sort(new TrancheComparator<>());
        return tranches;
    }

    public static List<TruthSensitivityTranche> mergeAndConvertTranches(final TreeMap<Double, List<VQSLODTranche>> scatteredTranches, final List<Double> tsLevels, final VariantRecalibratorArgumentCollection.Mode mode) {
        List<VQSLODTranche> mergedTranches = new ArrayList<>();
        List<TruthSensitivityTranche> gatheredTranches = new ArrayList<>();

        //make a list of merged tranches of the same length
        for (final Double VQSLODlevel : scatteredTranches.descendingKeySet()) {
            mergedTranches.add(mergeAndConvertTranches(scatteredTranches.get(VQSLODlevel),mode));
        }

        //go through the list of merged tranches to select those that match most closely the tsLevels
        //don't assume tsLevels are sorted
        sort(tsLevels);
        //tranches should be sorted based on the file format and the fact that we read them into a list

        ListIterator<Double> tsIter= tsLevels.listIterator();
        double targetTS = tsIter.next();
        double sensitivityDelta = 100.0; //initialize to 100%
        double prevDelta;
        ListIterator<VQSLODTranche> trancheIter = mergedTranches.listIterator();
        VQSLODTranche currentTranche = trancheIter.next();
        VQSLODTranche prevTranche;

        //match the calculated tranches with the requested truth sensitivity outputs as best as possible
        while (trancheIter.hasNext()) {  //we can only add tranches that we have and we should only add each at most once
            prevDelta = sensitivityDelta;
            prevTranche = currentTranche;
            currentTranche = trancheIter.next();
            sensitivityDelta = abs(targetTS - currentTranche.getTruthSensitivity()*100);
            if (sensitivityDelta > prevDelta) {
                gatheredTranches.add(new TruthSensitivityTranche(targetTS, prevTranche.minVQSLod, prevTranche.numKnown, prevTranche.knownTiTv, prevTranche.numNovel, prevTranche.novelTiTv, prevTranche.accessibleTruthSites, prevTranche.callsAtTruthSites, mode));
                if (tsIter.hasNext()) {
                    targetTS = tsIter.next();
                    sensitivityDelta = abs(targetTS - currentTranche.getTruthSensitivity()*100);
                }
                else {
                    break;
                }
            }
            if ( !trancheIter.hasNext()) { //if we haven't seen the best match, but we ran out of tranches
                gatheredTranches.add(new TruthSensitivityTranche(targetTS, currentTranche.minVQSLod, currentTranche.numKnown, currentTranche.knownTiTv, currentTranche.numNovel, currentTranche.novelTiTv, currentTranche.accessibleTruthSites, currentTranche.callsAtTruthSites, mode));
            }
        }

        return gatheredTranches;
    }

    public static VQSLODTranche mergeAndConvertTranches(final List<VQSLODTranche> scatteredTranches, VariantRecalibratorArgumentCollection.Mode mode) {
        double indexVQSLOD = scatteredTranches.get(0).minVQSLod;
        int sumNumKnown = 0;
        double sumKnownTransitions = 0;
        double sumKnownTransversions = 0;
        int sumNumNovel = 0;
        double sumNovelTransitions = 0;
        double sumNovelTransversions = 0;
        int sumAccessibleTruthSites = 0;
        int sumCallsAtTruthSites = 0;

        for (final VQSLODTranche tranche : scatteredTranches) {
            if (tranche.minVQSLod != indexVQSLOD)
                throw new IllegalStateException("Scattered tranches do not contain the same VQSLODs");
            sumNumKnown += tranche.numKnown;
            double trancheKnownTransitions = (tranche.knownTiTv*tranche.numKnown) / (1+tranche.knownTiTv);
            sumKnownTransitions += trancheKnownTransitions;
            sumKnownTransversions += (tranche.numKnown - trancheKnownTransitions);
            sumNumNovel += tranche.numNovel;
            double trancheNovelTransitions = (tranche.novelTiTv*tranche.numNovel) / (1+tranche.novelTiTv);
            sumNovelTransitions += trancheNovelTransitions;
            sumNovelTransversions += (tranche.numNovel - trancheNovelTransitions);
            sumAccessibleTruthSites += tranche.accessibleTruthSites;
            sumCallsAtTruthSites += tranche.callsAtTruthSites;
        }

        return new VQSLODTranche(indexVQSLOD, sumNumKnown, sumKnownTransitions/sumKnownTransversions, sumNumNovel, sumNovelTransitions/sumNovelTransversions, sumAccessibleTruthSites, sumCallsAtTruthSites, mode, "gathered" + indexVQSLOD);
    }
}
