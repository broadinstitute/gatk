package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.ClippingTail;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;


/**
 * Rank Sum Test of REF versus ALT base quality scores
 *
 * <p>This variant-level annotation tests compares the base qualities of the data supporting the reference allele with those supporting the alternate allele. The ideal result is a value close to zero, which indicates there is little to no difference. A negative value indicates that the bases supporting the alternate allele have lower quality scores than those supporting the reference allele. Conversely, a positive value indicates that the bases supporting the alternate allele have higher quality scores than those supporting the reference allele. Finding a statistically significant difference either way suggests that the sequencing process may have been biased or affected by an artifact.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test for base qualities (bases supporting REF vs. bases supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The base quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Rank sum test of REF versus ALT base quality scores (BaseQRankSum)")
public final class BaseQualityRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.BASE_QUAL_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc) {
        return getReadBaseQuality(read, refLoc);
    }

    public static OptionalDouble getReadBaseQuality(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);

        final int readCoordinate = ReadUtils.getReadCoordinateForReferenceCoordinateUpToEndOfRead(read, refLoc, ClippingTail.RIGHT_TAIL, true);
        return readCoordinate == ReadUtils.CLIPPING_GOAL_NOT_REACHED || readCoordinate < 0 || readCoordinate >= read.getLength() ? OptionalDouble.empty() : OptionalDouble.of(read.getBaseQuality(readCoordinate));
    }
}
