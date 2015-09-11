package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;

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
public final class Coverage extends InfoFieldAnnotation implements StandardAnnotation {

    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        Utils.nonNull(vc);
        if (perReadAlleleLikelihoodMap == null || perReadAlleleLikelihoodMap.isEmpty()) {
            return null;
        }

        final int depth = perReadAlleleLikelihoodMap.values().stream().mapToInt(maps -> maps.getLikelihoodReadMap().size()).sum();
        return Collections.singletonMap(getKeyNames().get(0), String.format("%d", depth));
    }

    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.DEPTH_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}
