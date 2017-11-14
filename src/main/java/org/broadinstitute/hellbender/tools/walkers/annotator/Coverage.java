package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;

import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Total depth of coverage per sample and over all samples.
 *
 * <p>This annotation is used to provide counts of read depth at two different levels, with some important differences. At the sample level (FORMAT), the DP value is the count of reads that passed the caller's internal quality control metrics (such as MAPQ > 17, for example). At the site level (INFO), the DP value is the unfiltered depth over all samples.</p>
 *
 * <p>See the method documentation on <a href="http://www.broadinstitute.org/gatk/guide/article?id=4721">using coverage information</a> for important interpretation details.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>If downsampling is enabled (as is done by default for some analyses to remove excessive coverage), the depth of coverage effectively seen by the caller may be inferior to the actual depth of coverage in the original file. If using "-dcov D", the maximum depth that can be seen for N samples will be N * D.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample</a></b> calculates depth of coverage for each allele per sample (AD).</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerSampleHC.php">DepthPerSampleHC</a></b> calculates depth of coverage after filtering by HaplotypeCaller.</li>
 * </ul>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Total depth of coverage per sample and over all samples (DP)")
public final class Coverage extends InfoFieldAnnotation implements StandardAnnotation, StandardMutectAnnotation {

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if (likelihoods == null || likelihoods.readCount() == 0) {
            return Collections.emptyMap();
        }

        final int depth = likelihoods.readCount();
        return Collections.singletonMap(getKeyNames().get(0), String.format("%d", depth));
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.DEPTH_KEY); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}
