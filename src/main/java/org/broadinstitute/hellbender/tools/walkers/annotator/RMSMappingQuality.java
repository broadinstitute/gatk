package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


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
public final class RMSMappingQuality extends InfoFieldAnnotation implements StandardAnnotation, ReducibleAnnotation {

    @Override
    public String getRawKeyName() { return GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY;}

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     */
    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap){
        Utils.nonNull(vc);
        if (perReadAlleleLikelihoodMap == null || perReadAlleleLikelihoodMap.isEmpty() ) {
            return null;
        }

        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<Number> myData = new ReducibleAnnotationData<>(null);
        calculateRawData(vc, perReadAlleleLikelihoodMap, myData);
        final String annotationString = formatedValue((double) myData.getAttributeMap().get(Allele.NO_CALL));
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public void calculateRawData(final VariantContext vc,
                                 final Map<String, PerReadAlleleLikelihoodMap> pralm,
                                 final ReducibleAnnotationData rawAnnotations){
        if (pralm.isEmpty()) {
            return;
        }

        //put this as a double, like GATK3.5
        final double squareSum = pralm.values().stream()
                .flatMap(likelihoodMap -> likelihoodMap.getReads().stream())
                .map(GATKRead::getMappingQuality)
                .filter(mq -> mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE).mapToDouble(mq -> mq * mq).sum();
        rawAnnotations.putAttribute(Allele.NO_CALL, squareSum);
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        return annotateRawData(ref, vc, perReadAlleleLikelihoodMap);
    }

    @VisibleForTesting
    static String formatedValue(double rms) {
        return String.format("%.2f", rms);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(getRawKeyName()); }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }
}