package org.broadinstitute.hellbender.utils.pileup;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.AlignmentAndReferenceContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PileupDetectionArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;


/**
 * Helper class for handling pileup allele detection supplement for assembly. This code is analogous but not exactly
 * equivalent to the DRAGEN ColumnwiseDetection approach.
 */
public final class PileupBasedAlleles {

    final static String MISMATCH_BASES_PERCENTAGE_TAG = "MZ";
    public static final double MISMATCH_BASES_PERCENTAGE_ADJUSMTENT = 1000.0;

    // internal values to be used for tagging the VCF info fields appropriately
    public static String PILEUP_ALLELE_SUPPORTING_READS = "PileupSupportingReads";
    public static String PILEUP_ALLELE_BAD_READS_TAG = "PileupSupportingBadReads";
    public static String PILEUP_ALLELE_ASSEMBLY_BAD_READS_TAG = "PileupAssemblyBadReads";
    public static String PILEUP_ALLELE_TOTAL_READS = "PileupAssemblyTotalReads";

    private final static int COUNT_IDX = 0;
    private final static int BAD_COUNT_IDX = 1;
    private final static int ASSEMBLY_BAD_COUNT_IDX = 2;

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
    public static List<VariantContext> getPileupVariantContexts(final List<AlignmentAndReferenceContext> alignmentAndReferenceContextList, final PileupDetectionArgumentCollection args, final SAMFileHeader headerForReads, final int minBaseQualityScore) {

        final List<VariantContext> pileupVariantList = new ArrayList<>();

        // Iterate over every base
        for (int i = 0; i < alignmentAndReferenceContextList.size(); i++) {
            AlignmentAndReferenceContext alignmentAndReferenceContext = alignmentAndReferenceContextList.get(i);
            boolean onlyTrackDeletions = false;
            //Skip all work on sites that aren't active according to our heuristic
            if (args.activeRegionPhredThreshold > 0.0 && alignmentAndReferenceContext.getActivityScore() < args.activeRegionPhredThreshold ) {
                //This solves a discordance with Illumina where they count deletions (and thus construct them and filter on the threshold) not at the anchor base but at the first
                //deleted base on the reference. Consequently we must allow deletions to be detected one base upstream of adctive regions.
                if (args.detectIndels && i+1 < alignmentAndReferenceContextList.size() && alignmentAndReferenceContextList.get(i+1).getActivityScore() > args.activeRegionPhredThreshold ) {
                    onlyTrackDeletions = true;
                } else {
                    continue;
                }
            }

            final AlignmentContext alignmentContext = alignmentAndReferenceContext.getAlignmentContext();
            final ReferenceContext referenceContext = alignmentAndReferenceContext.getReferenceContext();
            int numOfBases = alignmentContext.size();
            final ReadPileup pileup = alignmentContext.getBasePileup();
            final byte refBase = referenceContext.getBase();

            //Key for counts arrays [support, bad reads, assembly bad reads]
            Map<String, int[]> insertionCounts = new HashMap<>();
            Map<Integer, int[]> deletionCounts = new HashMap<>();
            Map<Byte, int[]> altCounts = new HashMap<>();

            for (PileupElement element : pileup) {
                final byte eachBase = element.getBase();

                // Subtract out low quality bases to mimic the reading active region determination //TODO this might need to also ignore the qual basees
                if (element.getQual() < minBaseQualityScore) {
                    numOfBases--;
                }

                // check to see that the base is not ref (and non-deletion) and increment the alt counts (and evaluate if the read is "bad")
                if (refBase != eachBase && eachBase != 'D' && element.getQual() > args.qualityForSnpsInPileupDetection) {
                    incrementAltCount(eachBase, altCounts,
                            evaluateBadRead(element.getRead(), referenceContext, args, headerForReads),
                            evaluateBadReadForAssembly(element.getRead(), referenceContext, args, headerForReads));
                }

                //NOTE: we count indels
                if (args.detectIndels) {
                    // now look for indels
                    if (element.isBeforeInsertion()) {
                        incrementInsertionCount(element.getBasesOfImmediatelyFollowingInsertion(), insertionCounts,
                                evaluateBadRead(element.getRead(), referenceContext, args, headerForReads),
                                evaluateBadReadForAssembly(element.getRead(), referenceContext, args, headerForReads));
                    }
                    if (element.isBeforeDeletionStart()) {
                        incrementDeletionCount(element.getLengthOfImmediatelyFollowingIndel(), deletionCounts,
                                evaluateBadRead(element.getRead(), referenceContext, args, headerForReads),
                                evaluateBadReadForAssembly(element.getRead(), referenceContext, args, headerForReads));

                    }
                }
            }

            // Evaluate the detected SNP alleles for this site
            if (!onlyTrackDeletions) {
                for (Map.Entry<Byte, int[]> allele : altCounts.entrySet()) {
                    List<Allele> alleles = new ArrayList<>();
                    alleles.add(Allele.create(referenceContext.getBase(), true));
                    alleles.add(Allele.create(allele.getKey()));
                    final VariantContextBuilder pileupSNP = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), alleles);
                    pileupVariantList.add(pileupSNP
                            .attribute(PILEUP_ALLELE_SUPPORTING_READS, allele.getValue()[COUNT_IDX])
                            .attribute(PILEUP_ALLELE_BAD_READS_TAG, allele.getValue()[BAD_COUNT_IDX])
                            .attribute(PILEUP_ALLELE_ASSEMBLY_BAD_READS_TAG, allele.getValue()[ASSEMBLY_BAD_COUNT_IDX])
                            .attribute(PILEUP_ALLELE_TOTAL_READS, numOfBases).make());
                }
            }

            if (args.detectIndels) {
                // Evaluate the detected Insertions alleles for this site
                if (!onlyTrackDeletions) {
                    for (Map.Entry<String, int[]> allele : insertionCounts.entrySet()) {
                        List<Allele> delAlleles = new ArrayList<>();
                        delAlleles.add(Allele.create(referenceContext.getBase(), true));
                        delAlleles.add(Allele.create((char) referenceContext.getBase() + allele.getKey()));
                        final VariantContextBuilder pileupInsertion = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd(), delAlleles);
                        pileupVariantList.add(pileupInsertion
                                .attribute(PILEUP_ALLELE_SUPPORTING_READS, allele.getValue()[COUNT_IDX])
                                .attribute(PILEUP_ALLELE_BAD_READS_TAG, allele.getValue()[BAD_COUNT_IDX])
                                .attribute(PILEUP_ALLELE_ASSEMBLY_BAD_READS_TAG, allele.getValue()[ASSEMBLY_BAD_COUNT_IDX])
                                .attribute(PILEUP_ALLELE_TOTAL_READS, numOfBases).make());
                    }
                }

                // Evaluate the detected Deletions alleles for this site
                for (Map.Entry<Integer, int[]> allele : deletionCounts.entrySet()) {
                    List<Allele> insAlleles = new ArrayList<>();
                    insAlleles.add(Allele.create(referenceContext.getBase(), false));
                    insAlleles.add(Allele.create(referenceContext.getBases(
                            new SimpleInterval(referenceContext.getContig(),
                                    alignmentContext.getStart(),
                                    alignmentContext.getEnd() + allele.getKey())),
                            true));
                    final VariantContextBuilder pileupInsertion = new VariantContextBuilder("pileup", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd() + allele.getKey(), insAlleles);
                    pileupVariantList.add(pileupInsertion
                            .attribute(PILEUP_ALLELE_SUPPORTING_READS, allele.getValue()[COUNT_IDX])
                            .attribute(PILEUP_ALLELE_BAD_READS_TAG, allele.getValue()[BAD_COUNT_IDX])
                            .attribute(PILEUP_ALLELE_ASSEMBLY_BAD_READS_TAG, allele.getValue()[ASSEMBLY_BAD_COUNT_IDX])
                            .attribute(PILEUP_ALLELE_TOTAL_READS, numOfBases).make());
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
    public static boolean passesFilters(final PileupDetectionArgumentCollection args, final VariantContext pileupVariant) {
        //Validation that this VC is the correct object type
        validatePileupVariant(pileupVariant);

        final int pileupSupport = pileupVariant.getAttributeAsInt(PILEUP_ALLELE_SUPPORTING_READS, 0);
        final int totalDepth = pileupVariant.getAttributeAsInt(PILEUP_ALLELE_TOTAL_READS, 0);
        return ((double) pileupSupport / (double) totalDepth) > (pileupVariant.isIndel() ? args.indelThreshold : args.snpThreshold)
                && totalDepth >= args.pileupAbsoluteDepth
                && ((args.badReadThreshold <= 0.0) || (double) pileupVariant.getAttributeAsInt(PILEUP_ALLELE_BAD_READS_TAG,0) / (double)pileupSupport <= args.badReadThreshold);
    }

    // TODO this is the most sketchy one... does a variant that fails pileup calling with only one bad read as support count as garbage by this tool...
    public static boolean shouldFilterAssemblyVariant(final PileupDetectionArgumentCollection args, final VariantContext pileupVariant) {
        //Validation that this VC is the correct object type
        validatePileupVariant(pileupVariant);
        final int pileupSupport = pileupVariant.getAttributeAsInt(PILEUP_ALLELE_SUPPORTING_READS, 0);
        final int assemblyBadReads = pileupVariant.getAttributeAsInt(PILEUP_ALLELE_ASSEMBLY_BAD_READS_TAG, 0);
        return ((args.assemblyBadReadThreshold > 0.0) && (double) assemblyBadReads / (double) pileupSupport >= args.assemblyBadReadThreshold);
    }

    private static void validatePileupVariant(final VariantContext pileupVariant) {
        Utils.nonNull(pileupVariant);
        if (pileupVariant.getAlleles().size() != 2 ||
                !pileupVariant.hasAttribute(PILEUP_ALLELE_ASSEMBLY_BAD_READS_TAG) ||
                !pileupVariant.hasAttribute(PILEUP_ALLELE_SUPPORTING_READS) ||
                !pileupVariant.hasAttribute(PILEUP_ALLELE_TOTAL_READS) ||
                !pileupVariant.hasAttribute(PILEUP_ALLELE_BAD_READS_TAG)) {
            throw new GATKException.ShouldNeverReachHereException("The supplied Variant Context "+pileupVariant.toString()+" is not a PileupVariantContext");
        }
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
        // Assert that the edit distance for the read is in line
        final Integer mismatchPercentage = read.getAttributeAsInteger(MISMATCH_BASES_PERCENTAGE_TAG);
        Utils.nonNull(mismatchPercentage);
        if ((mismatchPercentage / MISMATCH_BASES_PERCENTAGE_ADJUSMTENT) > args.badReadEditDistance) {
            return true;
        }

        //TODO add threshold descibed by illumina about insert size compared to the average
        if (args.templateLengthStd > 0 && args.templateLengthMean > 0) {
            SAMRecord samRecordForRead = read.convertToSAMRecord(headerForRead);
            int templateLength = samRecordForRead.getInferredInsertSize();
            // This is an illumina magic number... Its possible none of this is particularly important for Functional Equivalency.
            if (templateLength < args.templateLengthMean - 2.25 * args.templateLengthStd
                    || templateLength > args.templateLengthMean + 2.25 * args.templateLengthStd) {
                return true;
            }
        }
        return false;
    }

    @VisibleForTesting
    static boolean evaluateBadReadForAssembly(final GATKRead read, final ReferenceContext referenceContext, final PileupDetectionArgumentCollection args, final SAMFileHeader headerForRead) {
        if (args.assemblyBadReadThreshold <= 0.0) {
            return false;
        }
        // TODO other checks?
        Utils.nonNull(read.getAttributeAsInteger(MISMATCH_BASES_PERCENTAGE_TAG));
        return (read.getAttributeAsInteger(MISMATCH_BASES_PERCENTAGE_TAG) / MISMATCH_BASES_PERCENTAGE_ADJUSMTENT) > args.assemblyBadReadEditDistance;
    }

    // Helper methods to manage the badness counting arrays
    private static void incrementInsertionCount(String insertion, Map<String, int[]> insertionCounts, boolean bad, boolean assemblybad){
        int[] values = insertionCounts.computeIfAbsent(insertion, (i) -> new int[3] );
        values[COUNT_IDX]+=1; values[BAD_COUNT_IDX]+=bad?1:0; values[ASSEMBLY_BAD_COUNT_IDX]+=assemblybad?1:0;
    }
    private static void incrementDeletionCount(Integer deletion, Map<Integer, int[]> deletionCounts, boolean bad, boolean assemblybad){
        int[] values =  deletionCounts.computeIfAbsent(deletion, (i) -> new int[3] );
        values[COUNT_IDX]+=1; values[BAD_COUNT_IDX]+=bad?1:0; values[ASSEMBLY_BAD_COUNT_IDX]+=assemblybad?1:0;
    }
    private static void incrementAltCount(byte base, Map<Byte, int[]> altCounts, boolean bad, boolean assemblybad){
        int[] values = altCounts.computeIfAbsent(base, (i) -> new int[3] );
        values[COUNT_IDX]+=1; values[BAD_COUNT_IDX]+=bad?1:0; values[ASSEMBLY_BAD_COUNT_IDX]+=assemblybad?1:0;
    }


    public static void addMismatchPercentageToRead(final GATKRead read, final SAMFileHeader headerForRead, final ReferenceContext referenceContext) {
        //TODO this conversion is really unnecessary. Perhaps we should expose a new SequenceUtil like NM tag calculation?...
        if (read.hasAttribute(MISMATCH_BASES_PERCENTAGE_TAG)){
            return;
        }

        SAMRecord samRecordForRead = read.convertToSAMRecord(headerForRead);
        final int nmScore;
        if (! read.hasAttribute("NM")) {
            nmScore = SequenceUtil.calculateSamNmTag(samRecordForRead, referenceContext.getBases(new SimpleInterval(read)), read.getStart() - 1);
        } else {
            nmScore = read.getAttributeAsInteger("NM");
        }
        // We adjust the NM score by any indels in the read
        int adjustedNMScore = nmScore - read.getCigarElements().stream().filter(element -> element.getOperator().isIndel()).mapToInt(CigarElement::getLength).sum();

        // NOTE: we store the percentage as an integer x1000 because that is what the read attributes support. 4 units of precision should be more than enough for this ratio in the first place whereas 3 is probably not (reads over 100 bases get binned).
        read.setAttribute(MISMATCH_BASES_PERCENTAGE_TAG, (int) (MISMATCH_BASES_PERCENTAGE_ADJUSMTENT * adjustedNMScore / (read.getCigarElements().stream().filter(element -> element.getOperator().isAlignment()).mapToInt(CigarElement::getLength).sum() )));
    }
}
