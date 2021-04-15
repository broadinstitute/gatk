package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

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
                if (refBase != eachBase) {
                    incrementAltCount(eachBase, altCounts);
                }
            }

            List<Allele> alleles = new ArrayList<>();
            alleles.add(Allele.create(referenceContext.getBase(), true));
            // TODO: AH & BG add an option to deal with multiple alt alleles
            Optional<Map.Entry<Byte, Integer>> maxAlt = altCounts.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue));

            if (maxAlt.isPresent() && maxAlt.get().getValue() / numOfBases > 0.1) {
                alleles.add(Allele.create(maxAlt.get().getKey()));
                VariantContextBuilder pileupSNP = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), alleles);
                pileupSNPsList.add(pileupSNP.make());

            }
        }

        return pileupSNPsList;


    }

    private static void incrementAltCount(byte base, Map<Byte, Integer> altCounts){
        Byte baseObj = new Byte(base);
        if (!altCounts.containsKey(baseObj)){
            altCounts.put(baseObj,1);
        } else {
            altCounts.put(baseObj, altCounts.get(baseObj)+1);
        }
    }

}

/* Questions: how to get cigar info for each base?
How are insertions and deletions represented in Alignemnt context pileup.getbases
TODO: include Indels
* */