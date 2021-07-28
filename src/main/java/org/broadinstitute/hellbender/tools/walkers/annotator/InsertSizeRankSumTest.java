package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;

/**
 * Rank Sum Test for insert size differences between reads that support the REF and ALT alleles
 *
 * <p>This variant-level annotation tests whether there is evidence of bias in the fragment size of reads that support the reference and alternate alleles. </p>
 *
 * <p>A value close to zero indicates there is little to no difference in the fragment sizes of reads supporting the two alleles. A negative value indicates that reads that support the alternate allele tend to be shorter than reads that support the reference allele. Conversely, a positive value indicates reads that support the alternate allele tend to be longer than reads that support the reference allele. </p>
 *
 * <p>This annotation can be used ....</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test for the insert size of reads supporting REF vs. reads supporting ALT. See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The fragment size rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Rank sum test for fragments size of reads supporting REF versus ALT alleles (InsertSizeRankSum)")
public final class InsertSizeRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.INSERT_SIZE_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GATKRead read, final VariantContext vc) {
        return getInsertSize(read, vc);
    }

    public static OptionalDouble getInsertSize(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);
        return OptionalDouble.of(read.getFragmentLength());
    }
}
