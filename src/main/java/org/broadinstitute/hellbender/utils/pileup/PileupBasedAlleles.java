package org.broadinstitute.hellbender.utils.pileup;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.AlignmentAndReferenceContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PileupDetectionArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;


/**
 * Helper class for handling pileup allele detection supplement for assembly. This code is analogous but not exactly
 * equivalent to the DRAGEN ColumnwiseDetection approach.
 */
public final class PileupBasedAlleles {

    /**
     * Accepts the raw per-base pileups stored from the active region detection code and parses them for potential variants
     * that are visible in the pileups but might be dropped from assembly for any number of reasons. The basic algorithm works
     * as follows:
     *  - iterate over every pileup and count alt bases
     *      - (beta) detect insertions overlapping this site (CURRENTLY ONLY WORKS FOR INSERTIONS)
     *  - count "bad" reads as defined by Illumina filtering for pileup detection of variants {@Link #evaluateBadRead}
     *  - For each detected alt, evaluate if the number of alternate bases are sufficient to make the call and make a VariantContext.
     *
     * @param alignmentAndReferenceContextList  List of stored pileups and reference context information where every element is a base from the active region.
     *                                          NOTE: the expectation is that the stored pileups are based off of the ORIGINAL (un-clipped) reads from active region determination.
     * @param args                              Configuration arguments to use for filtering/annotations
     * @param headerForReads                    Header for the reads (only necessary for SAM file conversion)
     * @return A list of variant context objects corresponding to potential variants that pass our heuristics.
     */
    public static ArrayList<VariantContext> getPileupVariantContexts(final List<AlignmentAndReferenceContext> alignmentAndReferenceContextList, final PileupDetectionArgumentCollection args, final SAMFileHeader headerForReads) {

        final ArrayList<VariantContext> pileupVariantList = new ArrayList<>();

        // Iterate over every base
        for(AlignmentAndReferenceContext alignmentAndReferenceContext : alignmentAndReferenceContextList) {
            final AlignmentContext alignmentContext = alignmentAndReferenceContext.getAlignmentContext();
            final ReferenceContext referenceContext = alignmentAndReferenceContext.getReferenceContext();
            final int numOfBases = alignmentContext.size();
            final ReadPileup pileup = alignmentContext.getBasePileup();
            final byte refBase = referenceContext.getBase();

            Map<String, Integer> insertionCounts = new HashMap<>();

            Map<Byte, Integer> altCounts = new HashMap<>();

            int totalAltReads = 0;
            int totalAltBadReads = 0;

            for (PileupElement element : pileup) {
                final byte eachBase = element.getBase();

                // check to see that the base is not ref (and non-deletion) and increment the alt counts (and evaluate if the read is "bad")
                if (refBase != eachBase && eachBase != 'D') {
                    incrementAltCount(eachBase, altCounts);
                    totalAltReads++;
                    // Handle the "badness"
                    if (evaluateBadRead(element.getRead(), referenceContext, args, headerForReads)) {
                        totalAltBadReads++;
                    }
                }

                // TODO currently this only handles Insertions.
                if (args.detectIndels) {
                    // now look for indels
                    if (element.isBeforeInsertion()) {
                        incrementInsertionCount(element.getBasesOfImmediatelyFollowingInsertion(), insertionCounts);
                    }

                    //TODO this is possibly double dipping if there are snps adjacent to indels?
                    totalAltReads++;
                    // Handle the "badness"
                    if (evaluateBadRead(element.getRead(), referenceContext, args, headerForReads)) {
                        totalAltBadReads++;
                    }
                }

            }

            // Evaluate the detected SNP alleles for this site
            final List<Allele> alleles = new ArrayList<>();
            alleles.add(Allele.create(referenceContext.getBase(), true));
            final Optional<Map.Entry<Byte, Integer>> maxAlt = altCounts.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue));
            if (maxAlt.isPresent()
                    && passesFilters(args, false, numOfBases, totalAltBadReads, totalAltReads, maxAlt.get())) {

                alleles.add(Allele.create(maxAlt.get().getKey()));
                final VariantContextBuilder pileupSNP = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), alleles);
                pileupVariantList.add(pileupSNP.make());
            }

            // Evaluate the detected INDEL alleles for this site
            if (args.detectIndels) {
                final List<Allele> indelAlleles = new ArrayList<>();
                indelAlleles.add(Allele.create(referenceContext.getBase(), true));
                final Optional<Map.Entry<String, Integer>> maxIns = insertionCounts.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue));
                if (maxIns.isPresent()
                        && passesFilters(args, true, numOfBases, totalAltBadReads, totalAltReads, maxIns.get())) {

                    indelAlleles.add(Allele.create((char)referenceContext.getBase() + maxIns.get().getKey()));
                    final VariantContextBuilder pileupInsertion = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), indelAlleles);
                    pileupVariantList.add(pileupInsertion.make());
                }
            }
        }

        return pileupVariantList;
    }

    /**
     * Apply the filters to discovered alleles
     * - Does it have greater than snpThreshold fraction of bases support in the pileups?
     * - Does it have greater than pileupAbsoluteDepth number of reads supporting it?
     * - Are the reads supporting alts at the site greater than badReadThreshold percent "good"? //TODO evaluate if this is worth doing on a per-allele basis or otherwise
     */
    private static boolean passesFilters(final PileupDetectionArgumentCollection args, boolean indel,  final int numOfBases, final int totalAltBadReads, final int totalAltReads, final Map.Entry<?, Integer> maxAlt) {
        return ((float) maxAlt.getValue() / (float) numOfBases) > (indel ? args.indelThreshold : args.snpThreshold)
                && numOfBases >= args.pileupAbsoluteDepth
                && ((args.badReadThreshold <= 0.0) || (float) totalAltBadReads / (float)totalAltReads <= args.badReadThreshold);
    }

    /**
     * Based on the illumina PileupDetection filtering code: We apply a number of configurable heuristics to the reads that support
     * alt alleles that may be added and for each read evaluate if its "bad" or not by each of the heurisitcs. Currently they are:
     * - Secondary/SA tag reads are bad
     * - Improperly paired reads are bad
     * - Reads with > 8% per-base edit distance to the reference are bad
     * - Reads 2 std deviations away from the standard insert size are bad (not implemented)
     *
     * @param read
     * @param referenceContext
     * @param args
     * @param headerForRead TODO get rid of this sam record conversion
     * @return true if any of the "badness" heuristics suggest we should consider the read suspect, false otherwise.
     */
    @VisibleForTesting
    static boolean evaluateBadRead(final GATKRead read, final ReferenceContext referenceContext, final PileupDetectionArgumentCollection args, final SAMFileHeader headerForRead) {
        if (args.badReadThreshold <= 0.0) {
            return false;
        }
        if (args.badReadProperPair && !read.isProperlyPaired()) {
            return true;
        }
        if (args.badReadSecondaryOrSupplementary && (read.isSecondaryAlignment() || read.hasAttribute("SA"))) {
            return true;
        }


        //TODO this conversion is really unnecessary. Perhaps we should expose a new SequenceUtil like NM tag calculation?...
        SAMRecord samRecordForRead = read.convertToSAMRecord(headerForRead);

        // Assert that the edit distance for the read is in line
        if (args.badReadEditDistance > 0.0) {
            final int nmScore;
            if (! read.hasAttribute("NM")) {
                nmScore = SequenceUtil.calculateSamNmTag(samRecordForRead, referenceContext.getBases(new SimpleInterval(read)), read.getStart() - 1);
            } else {
                nmScore = read.getAttributeAsInteger("NM");
            }
            if (nmScore > (read.getLength() * args.badReadEditDistance)) {
                return true;
            }
        }

        //TODO add threshold descibed by illumina about insert size compared to the average
        if (args.templateLengthStd > 0 && args.templateLengthMean > 0) {
            int templateLength = samRecordForRead.getInferredInsertSize();
            // This is an illumina magic number... Its possible none of this is particularly important for Functional Equivalency.
            if (templateLength < args.templateLengthMean - 2.25 * args.templateLengthStd
                    || templateLength > args.templateLengthMean + 2.25 * args.templateLengthStd) {
                return true;
            }
        }
        return false;
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
