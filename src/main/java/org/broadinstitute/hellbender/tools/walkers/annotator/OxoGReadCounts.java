package org.broadinstitute.hellbender.tools.walkers.annotator;

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
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.utils.BaseUtils.Base.A;
import static org.broadinstitute.hellbender.utils.BaseUtils.Base.C;


/**
 *  Count of read pairs in the F1R2 and F2R1 configurations supporting the reference and alternate alleles
 *
 *  <p>This is an annotation that gathers information about the read pair configuration for the reads supporting each
 *  allele. It can be used along with downstream filtering steps to identify and filter out erroneous variants that occur
 *  with higher frequency in one read pair orientation.</p>
 *
 *  <h3>References</h3>
 *  <p>For more details about the mechanism of oxoG artifact generation, see <a href='http://www.ncbi.nlm.nih.gov/pubmed/23303777' target='_blank'>
 *      "Discovery and characterization of artefactual mutations in deep coverage targeted capture sequencing data due to oxidative DNA damage during sample preparation."
 *  by Costello et al.</a></p>
 */
public final class OxoGReadCounts extends GenotypeAnnotation {

    public static final Allele REF_C = Allele.create(C.base, true);
    public static final Allele REF_A = Allele.create(A.base, true);

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.OXOG_ALT_F1R2_KEY,
                             GATKVCFConstants.OXOG_ALT_F2R1_KEY,
                             GATKVCFConstants.OXOG_REF_F1R2_KEY,
                             GATKVCFConstants.OXOG_REF_F2R1_KEY,
                             GATKVCFConstants.OXOG_FRACTION_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_ALT_F1R2_KEY),
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_ALT_F2R1_KEY),
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_REF_F1R2_KEY),
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_REF_F2R1_KEY),
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.OXOG_FRACTION_KEY));
    }

    @Override
    public void annotate(final ReferenceContext refContext,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final ReadLikelihoods<Allele> likelihoods){
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        if (g == null || !g.isCalled() || likelihoods == null || !vc.isSNP()) {
            return;
        }

        final Allele ref = vc.getReference();
        final Allele alt = vc.getAlternateAllele(0);

        int alt_F1R2 = 0;
        int alt_F2R1 = 0;
        int ref_F1R2 = 0;
        int ref_F2R1 = 0;

        for (final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles(g.getSampleName())) {
            final GATKRead read = bestAllele.read;
            if (bestAllele.isInformative() && isUsableRead(read) && read.isPaired()) {
                final Allele allele = bestAllele.allele;
                if (allele.equals(ref, true)) {
                    if (read.isReverseStrand() == read.isFirstOfPair()) {
                        ref_F2R1++;
                    } else {
                        ref_F1R2++;
                    }
                } else if (allele.equals(alt, true)) {
                    if (read.isReverseStrand() == read.isFirstOfPair()) {
                        alt_F2R1++;
                    } else {
                        alt_F1R2++;
                    }
                }
            }
        }

        final double numerator;
        if (ref.equals(REF_C) || ref.equals(REF_A)) {
            numerator = alt_F2R1;
        } else {
            numerator = alt_F1R2;
        }
        final double denominator =  alt_F1R2 + alt_F2R1;
        final double fraction = numerator/denominator;

        gb.attribute(GATKVCFConstants.OXOG_ALT_F1R2_KEY, alt_F1R2);
        gb.attribute(GATKVCFConstants.OXOG_ALT_F2R1_KEY, alt_F2R1);
        gb.attribute(GATKVCFConstants.OXOG_REF_F1R2_KEY, ref_F1R2);
        gb.attribute(GATKVCFConstants.OXOG_REF_F2R1_KEY, ref_F2R1);
        gb.attribute(GATKVCFConstants.OXOG_FRACTION_KEY, fraction);
    }

    protected static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }
}