package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


/**
 * Root Mean Square of the mapping quality of reads across all samples.
 *
 * <p>This annotation provides an estimation of the overall mapping quality of reads supporting a variant call, averaged over all samples in a cohort.</p>
 *
 * <h3>Statistical notes</h3>
 * <p>The root mean square is equivalent to the mean of the mapping qualities plus the standard deviation of the mapping qualities.</p>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityRankSumTest.php">MappingQualityRankSumTest</a></b> compares the mapping quality of reads supporting the REF and ALT alleles.</li>
 * </ul>
 *
 */
public final class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation {

    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        Utils.nonNull(vc);
        if (perReadAlleleLikelihoodMap == null || perReadAlleleLikelihoodMap.isEmpty() ) {
            return null;
        }

        //get all available mapping qualities
        final List<Integer> qualities = perReadAlleleLikelihoodMap.values().stream().
                flatMap(pprl -> pprl.getReads().stream().map(read -> read.getMappingQuality()).filter(mq -> mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE))
                .collect(Collectors.toList());
        final double rms = MathUtils.rms(qualities);
        return Collections.singletonMap(getKeyNames().get(0), formatedValue(rms));
    }

    @VisibleForTesting
    static String formatedValue(double rms) {
        return String.format("%.2f", rms);
    }

    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.RMS_MAPPING_QUALITY_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}