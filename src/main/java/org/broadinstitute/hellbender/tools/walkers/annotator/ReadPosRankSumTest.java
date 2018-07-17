package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;

/**
 * Rank Sum Test for relative positioning of REF versus ALT alleles within reads
 *
 * <p>This variant-level annotation tests whether there is evidence of bias in the position of alleles within the reads that support them, between the reference and alternate alleles. Seeing an allele only near the ends of reads is indicative of error, because that is where sequencers tend to make the most errors. However, some variants located near the edges of sequenced regions will necessarily be covered by the ends of reads, so we can't just set an absolute "minimum distance from end of read" threshold. That is why we use a rank sum test to evaluate whether there is a difference in how well the reference allele and the alternate allele are supported.</p>
 *
 * <p>The ideal result is a value close to zero, which indicates there is little to no difference in where the alleles are found relative to the ends of reads. A negative value indicates that the alternate allele is found at the ends of reads more often than the reference allele. Conversely, a positive value indicates that the reference allele is found at the ends of reads more often than the alternate allele. </p>
 *
 * <p>This annotation can be used to evaluate confidence in a variant call and is a recommended covariate for variant recalibration (VQSR). Finding a statistically significant difference in relative position either way suggests that the sequencing process may have been biased or affected by an artifact. In practice, we only filter out low negative values when evaluating variant quality because the idea is to filter out variants for which the quality of the data supporting the alternate allele is comparatively low. The reverse case, where it is the quality of data supporting the reference allele that is lower (resulting in positive ranksum scores), is not really informative for filtering variants.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test for site position within reads (position within reads supporting REF vs. position within reads supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Rank sum test for relative positioning of REF versus ALT alleles within reads (ReadPosRankSum)")
public final class ReadPosRankSumTest extends RankSumTest implements StandardAnnotation {

    public static final int LEGACY_PAIRHMM_BASE_QUALITY_SCORE_THRESHOLD = 20; //TODO the 20 hard value here is in lieu of PairHMM.BASE_QUALITY_SCORE_THRESHOLD in order to directly match GATK3 output

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.READ_POS_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc) {
        return getReadPosition(read, refLoc);
    }

    @Override
    public boolean isUsableRead(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        return super.isUsableRead(read, refLoc) && read.getSoftEnd() >= refLoc;
    }

    public static OptionalDouble getReadPosition(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            return OptionalDouble.empty();
        }

        // If the offset inside a deletion, it does not lie on a read.
        if ( AlignmentUtils.isInsideDeletion(read.getCigar(), offset) ) {
            return OptionalDouble.of(INVALID_ELEMENT_FROM_READ);
        }

        // hard clips at this point in the code are perfectly good bases that were clipped to make the read fit the assembly region
        final Cigar cigar = read.getCigar();
        final CigarElement firstElement = cigar.getFirstCigarElement();
        final CigarElement lastElement = cigar.getLastCigarElement();
        final int leadingHardClips = firstElement.getOperator() == CigarOperator.HARD_CLIP ? firstElement.getLength() : 0;
        final int trailingHardClips = lastElement.getOperator() == CigarOperator.HARD_CLIP ? lastElement.getLength() : 0;

        int readPos = leadingHardClips + AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), offset, false, 0, 0);
        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );
        final int numOriginalBases = numAlignedBases + leadingHardClips + trailingHardClips;

        //After the middle of the read, we compute the postion from the end of the read.
        if (readPos > numOriginalBases / 2) {
            readPos = numOriginalBases - (readPos + 1);
        }
        return OptionalDouble.of(readPos);
    }

    @Override
    protected OptionalDouble getElementForPileupElement(final PileupElement p, int refLoc) {
        final int offset = AlignmentUtils.calcAlignmentByteArrayOffset(p.getRead().getCigar(), p, 0, 0);
        return OptionalDouble.of(getFinalVariantReadPosition(p.getRead(), offset));
    }

    // Utility methods necessary for computing the rank sum using read pileups.
    // TODO these methods contain bugs ported from GATK3, see issue https://github.com/broadinstitute/gatk/issues/4450 for progress towards fixing them

    /**
     * Get the position of a variant within a read with respect to the closer end, accounting for hard clipped bases and low quality ends
     * Used by ReadPosRankSum annotations
     *
     * @param read  a read containing the variant
     * @param initialReadPosition   the position based on the modified, post-hard-clipped CIGAR
     * @return read position
     */
    public static int getFinalVariantReadPosition(final GATKRead read, final int initialReadPosition) {
        final int numAlignedBases = getNumAlignedBases(read);

        int readPos = initialReadPosition;
        //TODO: this doesn't work for the middle-right position if we index from zero
        if (initialReadPosition > numAlignedBases / 2) {
            readPos = numAlignedBases - (initialReadPosition + 1);
        }
        return readPos;
    }

    /**
     *
     * @param read  read to count bases upon
     * @return  the number of hard clipped and low qual bases at the read start (where start is the leftmost end w.r.t. the reference)
     */
    public static int getNumClippedBasesAtStart(final GATKRead read) {
        // check for hard clips (never consider these bases):
        final CigarElement first = read.getCigarElement(0);

        final byte[] unclippedReadBases = read.getBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        //TODO: this code may not even get used because HaplotypeCaller already hard clips low quality tails

        int numStartClippedBasesCounter = (first.getOperator() == CigarOperator.H)? first.getLength(): 0;

        for (int i = numStartClippedBasesCounter; i < unclippedReadBases.length; i++) {
            if (unclippedReadQuals[i] < LEGACY_PAIRHMM_BASE_QUALITY_SCORE_THRESHOLD) {
                numStartClippedBasesCounter++;
            } else {
                break;
            }
        }

        return numStartClippedBasesCounter;
    }

    /**
     *
     * @param read  a read containing the variant
     * @return  number of non-hard clipped, aligned bases (excluding low quality bases at either end)
     */
    //TODO: this is bizarre -- this code counts hard clips, but then subtracts them from the read length, which already doesn't count hard clips
    public static int getNumAlignedBases(final GATKRead read) {
        return read.getLength() - getNumClippedBasesAtStart(read) - getNumClippedBasesAtEnd(read);
    }

    /**
     *
     * @param read  a read containing the variant
     * @return  number of hard clipped and low qual bases at the read end (where end is right end w.r.t. the reference)
     */
    public static int getNumClippedBasesAtEnd(final GATKRead read) {
        // check for hard clips (never consider these bases):
        CigarElement last = read.getCigarElement(read.numCigarElements() - 1);

        final byte[] unclippedReadBases = read.getBases();
        final byte[] unclippedReadQuals = read.getBaseQualities();

        // Do a stricter base clipping than provided by CIGAR string, since this one may be too conservative,
        // and may leave a string of Q2 bases still hanging off the reads.
        //TODO: this code may not even get used because HaplotypeCaller already hard clips low quality tails

        int endClippedBasesCounter = last.getOperator() == CigarOperator.H ? last.getLength() : 0;

        for (int i = unclippedReadBases.length - endClippedBasesCounter - 1; i >= 0; i--) {
            if (unclippedReadQuals[i] < LEGACY_PAIRHMM_BASE_QUALITY_SCORE_THRESHOLD) {

                endClippedBasesCounter++;
            } else {
                break;
            }
        }
        return endClippedBasesCounter;
    }
}
