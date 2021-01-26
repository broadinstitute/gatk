package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.mutect.AlignmentData;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import scala.Char;

import java.util.*;


public final class PileupBasedAlleles {


    ArrayList<VariantContext> forcedPileupAlleles;


    public static ArrayList<VariantContext> getPileupVariantContexts(List<AlignmentData> alignmentDataList) {

        ArrayList<VariantContext> pileupSNPsList = new ArrayList<>();

        for(AlignmentData alignmentData : alignmentDataList) {
            AlignmentContext alignmentContext = alignmentData.getAlignmentContext();
            ReferenceContext referenceContext = alignmentData.getReferenceContext();
            final int numOfBases = alignmentContext.size();
            final ReadPileup pileup = alignmentContext.getBasePileup();
            final byte refBase = referenceContext.getBase();

            Map<String, Integer> insertionCounts = new HashMap<>();

            Map<Byte, Integer> altCounts = new HashMap<>();

            for (PileupElement element : pileup) {
//            for (byte eachBase : pileup.getBases()) {
                byte eachBase = element.getBase();
                // check to see that the base is not ref and that the alleles are one of these bases - ATCGN
//                if (refBase != eachBase && Allele.acceptableAlleleBases(new byte[eachBase])) {
                // TODO: AH & BG add a better check for acceptable alleles for eachBase
                if (refBase != eachBase && eachBase != 68) {
                    incrementAltCount(eachBase, altCounts);
                }

                // now look for indels
                if (element.isBeforeInsertion()) {
                   incrementInsertionCount(element.getBasesOfImmediatelyFollowingInsertion(), insertionCounts);
                }

                // TODO BG & AH check for indels at the end of the active region using isBeforeInsertion
            }

            List<Allele> alleles = new ArrayList<>();
            alleles.add(Allele.create(referenceContext.getBase(), true));
            // TODO: AH & BG add an option to deal with multiple alt alleles
            Optional<Map.Entry<Byte, Integer>> maxAlt = altCounts.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue));
            if (maxAlt.isPresent() && ((float)maxAlt.get().getValue() / (float)numOfBases) > 0.10 && numOfBases >= 5 ) {
//            if (maxAlt.isPresent() && ((float)maxAlt.get().getValue() / (float)numOfBases) > 0.1 && maxAlt.get().getValue() > 10 ) {
                alleles.add(Allele.create(maxAlt.get().getKey()));
                VariantContextBuilder pileupSNP = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), alleles);
                pileupSNPsList.add(pileupSNP.make());
            }

            // evaluation indels
            List<Allele> indelAlleles = new ArrayList<>();
            indelAlleles.add(Allele.create(referenceContext.getBase(), true));
            Optional<Map.Entry<String, Integer>> maxIns = insertionCounts.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue));
            if (maxIns.isPresent() && ((float)maxIns.get().getValue() / (float)numOfBases) > 0.50 && numOfBases >= 5 ) {
                indelAlleles.add(Allele.create(String.format("%c", referenceContext.getBase()) + maxIns.get().getKey()));
                VariantContextBuilder pileupInsertion = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), indelAlleles);
                pileupSNPsList.add(pileupInsertion.make());
            }

        }

        return pileupSNPsList;


    }

    private static void incrementInsertionCount(String insertion, Map<String, Integer> insertionCounts){
        if (!insertionCounts.containsKey(insertion)) {
            insertionCounts.put(insertion, 1);
        } else {
            insertionCounts.put(insertion, insertionCounts.get(insertion) + 1);
        }
    }

    private static void incrementAltCount(byte base, Map<Byte, Integer> altCounts){
        Byte baseObj = Byte.valueOf(base);
            if (!altCounts.containsKey(baseObj)) {
                altCounts.put(baseObj, 1);
            } else {
                altCounts.put(baseObj, altCounts.get(baseObj) + 1);
            }
    }

}

/* Questions: how to get cigar info for each base?
How are insertions and deletions represented in Alignment context pileup.getbases
Args to test - error-correction-log-odds; this is not turned on by default; may help with precision
TODO: include Indels
* */
