package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQualityRankSumTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;


/**
 * Allele-specific rank Sum Test of REF versus ALT base quality scores
 *
 * <p>This variant-level annotation compares the base qualities of the data supporting the reference allele with those supporting each alternate allele. To be clear, it does so separately for each alternate allele. </p>
 *
 * <p>The ideal result is a value close to zero, which indicates there is little to no difference. A negative value indicates that the bases supporting the alternate allele have lower quality scores than those supporting the reference allele. Conversely, a positive value indicates that the bases supporting the alternate allele have higher quality scores than those supporting the reference allele. Finding a statistically significant difference either way suggests that the sequencing process may have been biased or affected by an artifact.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The value output for this annotation is the u-based z-approximation from the Mann-Whitney-Wilcoxon Rank Sum Test for base qualities (bases supporting REF vs. bases supporting ALT). See the <a href="http://www.broadinstitute.org/gatk/guide/article?id=4732">method document on statistical tests</a> for a more detailed explanation of the ranksum test.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>Uninformative reads are not used in these calculations.</li>
 * <li>The base quality rank sum test cannot be calculated for sites without a mixture of reads showing both the reference and alternate alleles.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <ul>
 * <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_BaseQualityRankSumTest.php">BaseQualityRankSumTest</a></b> outputs a version of this annotation that includes all alternate alleles in a single calculation.</li>
 * </ul>
 *
 */
public class AS_BaseQualityRankSumTest extends AS_RankSumTest implements AS_StandardAnnotation {

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.AS_BASE_QUAL_RANK_SUM_KEY);
    }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.AS_RAW_BASE_QUAL_RANK_SUM_KEY;}

    /**
     * Get the element for the given read at the given reference position
     *
     * @param read     the read
     * @param refLoc   the reference position
     * @return a Double representing the element to be used in the rank sum test, or null if it should not be used
     */
    @Override
    protected OptionalDouble getElementForRead(final GATKRead read, final int refLoc) {
        return BaseQualityRankSumTest.getReadBaseQuality(read, refLoc);
    }

}