package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

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
public final class ReadPosRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.READ_POS_RANK_SUM_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc) {
        return getReadPosition(read, refLoc);
    }

    @Override
    public boolean isUsableRead(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        return super.isUsableRead(read, refLoc) && ReadUtils.getSoftEnd(read) >= refLoc;
    }

    public static OptionalDouble getReadPosition(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        final int offset = ReadUtils.getReadCoordinateForReferenceCoordinate(ReadUtils.getSoftStart(read), read.getCigar(), refLoc, ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( offset == ReadUtils.CLIPPING_GOAL_NOT_REACHED ) {
            return OptionalDouble.empty();
        }

        // If the offset inside a deletion, it does not lie on a read.
        if ( AlignmentUtils.isInsideDeletion(read.getCigar(), offset) ) {
            return OptionalDouble.of(INVALID_ELEMENT_FROM_READ);
        }

        int readPos = AlignmentUtils.calcAlignmentByteArrayOffset(read.getCigar(), offset, false, 0, 0);
        final int numAlignedBases = AlignmentUtils.getNumAlignedBasesCountingSoftClips( read );

        //After the middle of the read, we compute the postion from the end of the read.
        if (readPos > numAlignedBases / 2) {
            readPos = numAlignedBases - (readPos + 1);
        }
        return OptionalDouble.of(readPos);
    }


}
