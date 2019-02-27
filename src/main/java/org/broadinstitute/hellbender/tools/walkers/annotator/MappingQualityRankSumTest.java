package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;
import java.util.OptionalDouble;


/**
 * Rank Sum Test for mapping qualities of REF versus ALT reads
 *
 * <p>This variant-level annotation compares the mapping qualities of the reads supporting the reference allele with those supporting the alternate allele. The ideal result is a value close to zero, which indicates there is little to no difference. A negative value indicates that the reads supporting the alternate allele have lower mapping quality scores than those supporting the reference allele. Conversely, a positive value indicates that the reads supporting the alternate allele have higher mapping quality scores than those supporting the reference allele.</p>
 * <p>This annotation can be used to evaluate confidence in a variant call and is a recommended covariate for variant recalibration (VQSR). Finding a statistically significant difference in quality either way suggests that the sequencing and/or mapping process may have been biased or affected by an artifact. In practice, we only filter out low negative values when evaluating variant quality because the idea is to filter out variants for which the quality of the data supporting the alternate allele is comparatively low. The reverse case, where it is the quality of data supporting the reference allele that is lower (resulting in positive ranksum scores), is not really informative for filtering variants.
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test for mapping qualities (MAPQ of reads supporting REF vs. MAPQ of reads supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveat</h3>
 * <p>The mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_RMSMappingQuality.php">RMSMappingQuality</a></b> gives an estimation of the overal read mapping quality supporting a variant call.</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Rank sum test for mapping qualities of REF versus ALT reads (MQRankSum)")
public final class MappingQualityRankSumTest extends RankSumTest implements StandardAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY); }

    @Override
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc) {
        Utils.nonNull(read);
        return OptionalDouble.of(read.getMappingQuality());
    }

    protected OptionalDouble getElementForPileupElement(final PileupElement p, final int refLoc) {
        // default to returning the same value
        return OptionalDouble.of(p.getMappingQual());
    }
}
