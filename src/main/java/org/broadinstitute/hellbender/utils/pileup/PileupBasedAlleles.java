package org.broadinstitute.hellbender.utils.pileup;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Streams;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.AlignmentAndReferenceContext;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PileupDetectionArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;


/**
 * Helper class for handling pileup allele detection supplement for assembly. This code is analogous but not exactly
 * equivalent to the DRAGEN ColumnwiseDetection approach.
 */
public final class PileupBasedAlleles {

    final static String MISMATCH_BASES_PERCENTAGE_TAG = "MZ";
    public static final double MISMATCH_BASES_PERCENTAGE_ADJUSMTENT = 1000.0;

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
    public static Pair<Set<Event>, Set<Event>> goodAndBadPileupEvents(final List<AlignmentAndReferenceContext> alignmentAndReferenceContextList, final PileupDetectionArgumentCollection args, final SAMFileHeader headerForReads, final int minBaseQualityScore) {
        if (!args.usePileupDetection) {
            return ImmutablePair.of(Collections.emptySet(), Collections.emptySet());
        }

        final Set<Event> goodEvents = new HashSet<>();
        final Set<Event> badEvents = new HashSet<>();

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
            final String contig = alignmentContext.getContig();
            final int start = alignmentContext.getStart();
            final int end = alignmentContext.getEnd();

            final ReferenceContext referenceContext = alignmentAndReferenceContext.getReferenceContext();
            final MutableInt pileupDepth = new MutableInt(alignmentContext.size());
            final ReadPileup pileup = alignmentContext.getBasePileup();
            final byte refBase = referenceContext.getBase();

            //Key for counts arrays [support, bad reads, assembly bad reads]
            Map<String, int[]> insertionCounts = new HashMap<>();
            Map<Integer, int[]> deletionCounts = new HashMap<>();
            Map<Byte, int[]> SNPCounts = new HashMap<>();

            for (PileupElement element : pileup) {
                final byte eachBase = element.getBase();

                // Subtract out low quality bases to mimic the reading active region determination //TODO this might need to also ignore the qual basees
                if (element.getQual() < minBaseQualityScore) {
                    pileupDepth.decrement();
                }

                final boolean SNPFound = !onlyTrackDeletions && refBase != eachBase && eachBase != 'D' && element.getQual() > args.qualityForSnpsInPileupDetection;
                final boolean insertionFound = !onlyTrackDeletions && args.detectIndels && element.isBeforeInsertion();
                final boolean deletionFound = args.detectIndels && element.isBeforeDeletionStart();

                if (SNPFound || insertionFound || deletionFound) {
                    final boolean badPileup = badPileupRead(element.getRead(), args, headerForReads);
                    final boolean badAssembly = badAssemblyRead(element.getRead(), args);

                    if (SNPFound) {
                        incrementCounts(eachBase, SNPCounts, badPileup, badAssembly);
                    }

                    if (insertionFound) {
                        incrementCounts(element.getBasesOfImmediatelyFollowingInsertion(), insertionCounts, badPileup, badAssembly);
                    }

                    if (deletionFound) {
                        incrementCounts(element.getLengthOfImmediatelyFollowingIndel(), deletionCounts, badPileup, badAssembly);
                    }
                }
            }

            final Map<Event, int[]> SNPEventsAndCounts = SNPCounts.entrySet().stream()
                    .collect(Collectors.toMap(entry -> new Event(contig, start, Allele.create(refBase, true), Allele.create(entry.getKey())), entry -> entry.getValue()));
            final Map<Event, int[]> insertionEventsAndCounts = insertionCounts.entrySet().stream()
                    .collect(Collectors.toMap(entry -> new Event(contig, start, Allele.create(refBase, true), Allele.create((char) refBase + entry.getKey())), entry -> entry.getValue()));
            final Map<Event, int[]> deletionEventsAndCounts = deletionCounts.entrySet().stream()
                    .collect(Collectors.toMap(entry -> new Event(contig, start, Allele.create(referenceContext.getBases(new SimpleInterval(contig, start, end + entry.getKey())), true), Allele.create(refBase)), entry -> entry.getValue()));

            Streams.concat(SNPEventsAndCounts.entrySet().stream(), insertionEventsAndCounts.entrySet().stream(), deletionEventsAndCounts.entrySet().stream())
                    .forEach(eventAndCounts -> {
                        final Event event = eventAndCounts.getKey();
                        final int[] counts = eventAndCounts.getValue();

                        if (passesPileupFilters(args, counts[COUNT_IDX], counts[BAD_COUNT_IDX], pileupDepth.intValue(), event.isIndel())) {
                            goodEvents.add(event);
                        }

                        if (failsAssemblyFilters(args, counts[COUNT_IDX], counts[ASSEMBLY_BAD_COUNT_IDX])) {
                            badEvents.add(event);
                        }
                    });
        }

        return ImmutablePair.of(goodEvents, badEvents);
    }

    /**
     * Apply the filters to discovered alleles
     * - Does it have greater than snpThreshold fraction of bases support in the pileups?
     * - Does it have greater than pileupAbsoluteDepth number of reads supporting it?
     * - Are the reads supporting alts at the site greater than badReadThreshold percent "good"? //TODO evaluate if this is worth doing on a per-allele basis or otherwise
     */
    private static boolean passesPileupFilters(final PileupDetectionArgumentCollection args, final int pileupSupport, final int pileupBadReads, final int pileupDepth, final boolean isIndel) {
        return ((double) pileupSupport / (double) pileupDepth) > (isIndel ? args.indelThreshold : args.snpThreshold)
                && pileupDepth >= args.pileupAbsoluteDepth
                && ((args.badReadThreshold <= 0.0) || (double) pileupBadReads / (double)pileupSupport <= args.badReadThreshold);
    }

    // TODO this is the most sketchy one... does a variant that fails pileup calling with only one bad read as support count as garbage by this tool...
    private static boolean failsAssemblyFilters(final PileupDetectionArgumentCollection args, final int pileupSupport, final int assemblyBadReads) {
        return (args.assemblyBadReadThreshold > 0.0) && (double) assemblyBadReads / (double) pileupSupport >= args.assemblyBadReadThreshold;
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
     * @param args
     * @param headerForRead TODO get rid of this sam record conversion
     * @return true if any of the "badness" heuristics suggest we should consider the read suspect, false otherwise.
     */
    @VisibleForTesting
    static boolean badPileupRead(final GATKRead read, final PileupDetectionArgumentCollection args, final SAMFileHeader headerForRead) {
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
    static boolean badAssemblyRead(final GATKRead read, final PileupDetectionArgumentCollection args) {
        if (args.assemblyBadReadThreshold <= 0.0) {
            return false;
        }
        // TODO other checks?
        Utils.nonNull(read.getAttributeAsInteger(MISMATCH_BASES_PERCENTAGE_TAG));
        return (read.getAttributeAsInteger(MISMATCH_BASES_PERCENTAGE_TAG) / MISMATCH_BASES_PERCENTAGE_ADJUSMTENT) > args.assemblyBadReadEditDistance;
    }

    // Helper method to manage the badness counting arrays
    // T is a Byte for a SNP (alt base), String for insertion (inserted bases), Integer for deletion (deletion length)
    private static <T> void incrementCounts(T altAllele, Map<T, int[]> altCounts, boolean pileupBad, boolean assemblyBad){
        final int[] values = altCounts.computeIfAbsent(altAllele, (i) -> new int[3] );
        values[COUNT_IDX]+=1; values[BAD_COUNT_IDX]+=pileupBad?1:0; values[ASSEMBLY_BAD_COUNT_IDX]+=assemblyBad?1:0;
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
