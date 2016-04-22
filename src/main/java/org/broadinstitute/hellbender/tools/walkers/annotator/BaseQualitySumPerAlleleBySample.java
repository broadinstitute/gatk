package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;
import java.util.Map;


/**
 * Sum of evidence in reads supporting each allele for each sample
 *
 * <p>In the domain of somatic variants, a variant call can be supported by a few high quality reads. The
 * BaseQualitySumPerAlleleBySample annotation aims to give the user an estimate of the quality of the evidence supporting
 * a variant.</p>
 *
 * <h3>Notes</h3>
 * BaseQualitySumPerAlleleBySample is called and used by M2 for variant filtering. This annotation is applied to SNPs
 * and INDELs. Qualities are not literal base qualities, but instead are derived from the per-allele likelihoods derived
 * from the assembly engine.
 */
public final class BaseQualitySumPerAlleleBySample extends GenotypeAnnotation implements StandardSomaticAnnotation {

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(GATKVCFConstants.QUALITY_SCORE_SUM_KEY); }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getFormatLine(getKeyNames().get(0)));
    }

    @Override
    public void annotate(final ReferenceContext ref,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final PerReadAlleleLikelihoodMap alleleLikelihoodMap) {
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");
        if ( g == null || !g.isCalled() || alleleLikelihoodMap == null ) {
            return;
        }

        double refQualSum= 0.0;
        double altQualSum= 0.0;

        for ( final Map.Entry<GATKRead, Map<Allele,Double>> el : alleleLikelihoodMap.getLikelihoodReadMap().entrySet() ) {
            final MostLikelyAllele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if ( ! a.isInformative() ) {
                continue; // read is non-informative
            }

            final GATKRead read = el.getKey();
            if ( isUsableRead(read) ) {
                final double baseQual = getBaseQualityForRead(read, vc.getStart());
                if ( a.getMostLikelyAllele().isReference() ) {
                    refQualSum += baseQual;
                } else if ( vc.getAlleles().contains(a.getMostLikelyAllele()) ) {
                    altQualSum += baseQual;
                }
            }
        }

        gb.attribute(GATKVCFConstants.QUALITY_SCORE_SUM_KEY, new Integer[]{(int)Math.floor(refQualSum), (int) Math.floor(altQualSum)});
    }

    @VisibleForTesting
    static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    private static double getBaseQualityForRead(final GATKRead read, final int refLoc) {
        return read.getBaseQualities()[ReadUtils.getReadCoordinateForReferenceCoordinateUpToEndOfRead(read, refLoc, ReadUtils.ClippingTail.RIGHT_TAIL)];
    }
}