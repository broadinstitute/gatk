package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.StrandBiasUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;
import java.util.Map;

public class AS_StrandBiasMutectAnnotation extends InfoFieldAnnotation implements StandardMutectAnnotation {
    private final static Logger logger = LogManager.getLogger(StrandBiasBySample.class);
    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);

        if ( likelihoods == null ) {
            logger.warn("Annotation will not be calculated, alleleLikelihoodMap is null");
            return null;
        }

        return StrandBiasUtils.computeSBAnnotation(vc, likelihoods, GATKVCFConstants.AS_SB_TABLE_KEY);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        return super.getDescriptions();
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.AS_SB_TABLE_KEY);
    }
}
