package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;

/**
 * Rank Sum Test for hard-clipped bases on REF versus ALT reads
 *
 * <p>This variant-level annotation tests whether the data supporting the reference allele shows more or less base clipping (hard clips) than those supporting the alternate allele. The ideal result is a value close to zero, which indicates there is little to no difference.  A negative value indicates that the reads supporting the alternate allele have more hard-clipped bases than those supporting the reference allele. Conversely, a positive value indicates that the reads supporting the alternate allele have fewer hard-clipped bases than those supporting the reference allele. Finding a statistically significant difference either way suggests that the sequencing and/or mapping process may have been biased or affected by an artifact.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test applied to base clips (number of hard-clipped bases on reads supporting REF vs. number of hard-clipped bases on reads supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The clipping rank sum test cannot be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 * <h3>Really Big Caveat</h3>
 * <p> In AssemblyRegionWalkers the annotation engine receives reads after they have been hard-clipped to fit the assembly region.
 *  Thus this annotation should not be used with HaplotypeCaller and Mutect2.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Rank sum test for hard-clipped bases on REF versus ALT reads (ClippingRankSum)")
public final class ClippingRankSumTest extends RankSumTest {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.CLIPPING_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        return OptionalDouble.of(AlignmentUtils.getNumHardClippedBases(read));
    }

    @Override
    protected OptionalDouble getElementForPileupElement(final PileupElement p, final int refLoc) {
        Utils.nonNull(p);
        // default to returning the same value
        return getElementForRead(p.getRead(),refLoc);
    }
 }
