package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import scala.Char;

import java.util.*;


final class PileupBasedAlleles {


    ArrayList<VariantContext> forcedPileupAlleles;


    public static ArrayList<VariantContext> getPileupVariantContexts(List<AlignmentData> alignmentDataList) {

        ArrayList<VariantContext> pileupSNPsList = new ArrayList<>();

        for(AlignmentData alignmentData : alignmentDataList) {
            AlignmentContext alignmentContext = alignmentData.getAlignmentContext();
            ReferenceContext referenceContext = alignmentData.getReferenceContext();
            final int numOfBases = alignmentContext.size();
            final ReadPileup pileup = alignmentContext.getBasePileup();
            final byte refBase = referenceContext.getBase();

            Map<Byte, Integer> altCounts = new HashMap<>();

            for (byte eachBase : pileup.getBases()) {
                // check to see that the base is not ref and that the alleles are one of these bases - ATCGN
//                if (refBase != eachBase && Allele.acceptableAlleleBases(new byte[eachBase])) {
                // TODO: AH & BG add a better check for acceptable alleles for eachBase
                if (refBase != eachBase && eachBase != 68) {
                    incrementAltCount(eachBase, altCounts);
                }
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
        }

        return pileupSNPsList;


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