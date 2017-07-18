package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.apache.commons.lang.mutable.MutableInt;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.F1R2_KEY, GATKVCFConstants.F2R1_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F1R2_KEY),
                GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.F2R1_KEY));
    }

    @Override
    public void annotate(final ReferenceContext refContext,
                                  final VariantContext vc,
                                  final Genotype g,
                                  final GenotypeBuilder gb,
                                  final ReadLikelihoods<Allele> likelihoods){
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        if (g == null || likelihoods == null) {
            return;
        }

        final Map<Allele, MutableInt> f1r2Counts = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        final Map<Allele, MutableInt> f2r1Counts = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new MutableInt(0)));

        Utils.stream(likelihoods.bestAlleles(g.getSampleName()))
                .filter(ba -> ba.isInformative() && isUsableRead(ba.read))
                .forEach(ba -> (isF2R1(ba.read) ? f2r1Counts : f1r2Counts).get(ba.allele).increment());

        final int[] f1r2 = likelihoods.alleles().stream().mapToInt(a -> f1r2Counts.get(a).intValue()).toArray();

        final int[] f2r1 = likelihoods.alleles().stream().mapToInt(a -> f2r1Counts.get(a).intValue()).toArray();

        gb.attribute(GATKVCFConstants.F1R2_KEY, f1r2);
        gb.attribute(GATKVCFConstants.F2R1_KEY, f2r1);
    }

    protected static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    protected static boolean isF2R1(final GATKRead read) {
        return read.isReverseStrand() == read.isFirstOfPair();
    }
}