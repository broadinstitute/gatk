package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.AlignmentAndReferenceContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;


public final class PileupBasedAlleles {


//    ArrayList<VariantContext> forcedPileupAlleles;


    public static ArrayList<VariantContext> getPileupVariantContexts(List<AlignmentAndReferenceContext> alignmentAndReferenceContextList) {

        final ArrayList<VariantContext> pileupSNPsList = new ArrayList<>();

        for(AlignmentAndReferenceContext alignmentAndReferenceContext : alignmentAndReferenceContextList) {
            final AlignmentContext alignmentContext = alignmentAndReferenceContext.getAlignmentContext();
            final ReferenceContext referenceContext = alignmentAndReferenceContext.getReferenceContext();
            final int numOfBases = alignmentContext.size();
            final ReadPileup pileup = alignmentContext.getBasePileup();
            final byte refBase = referenceContext.getBase();

            Map<String, Integer> insertionCounts = new HashMap<>();

            Map<Byte, Integer> altCounts = new HashMap<>();

            for (PileupElement element : pileup) {
                final byte eachBase = element.getBase();

                // check to see that the base is not ref and that the alleles are one of these bases - ATCGN
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

            final List<Allele> alleles = new ArrayList<>();
            alleles.add(Allele.create(referenceContext.getBase(), true));
            // TODO: AH & BG add an option to deal with multiple alt alleles
            final Optional<Map.Entry<Byte, Integer>> maxAlt = altCounts.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue));
            if (maxAlt.isPresent() && ((float)maxAlt.get().getValue() / (float)numOfBases) > 0.10 && numOfBases >= 5 ) {
                alleles.add(Allele.create(maxAlt.get().getKey()));
                final VariantContextBuilder pileupSNP = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), alleles);
                pileupSNPsList.add(pileupSNP.make());
            }

            // evaluation indels
            final List<Allele> indelAlleles = new ArrayList<>();
            indelAlleles.add(Allele.create(referenceContext.getBase(), true));
            final Optional<Map.Entry<String, Integer>> maxIns = insertionCounts.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue));
            if (maxIns.isPresent() && ((float)maxIns.get().getValue() / (float)numOfBases) > 0.50 && numOfBases >= 5 ) {
                indelAlleles.add(Allele.create(referenceContext.getBase() + maxIns.get().getKey()));
                final VariantContextBuilder pileupInsertion = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), indelAlleles);
                pileupSNPsList.add(pileupInsertion.make());
            }
        }

        return pileupSNPsList;
    }

    private static void incrementInsertionCount(String insertion, Map<String, Integer> insertionCounts){
        insertionCounts.put(insertion,
                insertionCounts.getOrDefault(insertion,0) + 1);
    }

    private static void incrementAltCount(byte base, Map<Byte, Integer> altCounts){
        altCounts.put(base,
                altCounts.getOrDefault(base,0) + 1);
    }
}

/* Questions: how to get cigar info for each base?
How are insertions and deletions represented in Alignment context pileup.getbases
Args to test - error-correction-log-odds; this is not turned on by default; may help with precision
TODO: include Indels
* */
