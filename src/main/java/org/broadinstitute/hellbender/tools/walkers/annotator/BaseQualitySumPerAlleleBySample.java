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
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;


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
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");
        if ( g == null || !g.isCalled() || likelihoods == null ) {
            return;
        }

        double refQualSum= 0.0;
        double altQualSum= 0.0;

        for ( final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles(g.getSampleName()) ) {
            if ( bestAllele.isInformative() && isUsableRead(bestAllele.read)) {
                final double baseQual = getBaseQualityForRead(bestAllele.read, vc.getStart());
                if (bestAllele.allele.isReference()) {
                    refQualSum += baseQual;
                } else if (vc.getAlleles().contains(bestAllele.allele)) {
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