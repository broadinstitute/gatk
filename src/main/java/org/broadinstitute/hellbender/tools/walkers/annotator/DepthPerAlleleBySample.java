package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;

import java.util.*;

/**
 * Depth of coverage of each allele per sample
 *
 * <p>Also known as the allele depth, this annotation gives the unfiltered count of reads that support a given allele for an individual sample. The values in the field are ordered to match the order of alleles specified in the REF and ALT fields: REF, ALT1, ALT2 and so on if there are multiple ALT alleles.</p>
 *
 * <p>See the method documentation on <a href="http://www.broadinstitute.org/gatk/guide/article?id=4721">using coverage information</a> for important interpretation details.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>The AD calculation as performed by HaplotypeCaller may not yield exact results because only reads that statistically favor one allele over the other are counted. Due to this fact, the sum of AD may be different than the individual sample depth, especially when there are many non-informative reads.</li>
 *     <li>Because the AD includes reads and bases that were filtered by the caller (and in case of indels, is based on a statistical computation), it  should not be used to make assumptions about the genotype that it is associated with. Ultimately, the phred-scaled genotype likelihoods (PLs) are what determines the genotype calls.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_Coverage.php">Coverage</a></b> gives the filtered depth of coverage for each sample and the unfiltered depth across all samples.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_AlleleBalance.php">AlleleBalance</a></b> is a generalization of this annotation over all samples.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_AlleleBalanceBySample.php">AlleleBalanceBySample</a></b> calculates allele balance for each individual sample.</li>
 * </ul>
 */
public final class DepthPerAlleleBySample extends GenotypeAnnotation implements StandardAnnotation {

    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final PerReadAlleleLikelihoodMap alleleLikelihoodMap) {
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        if ( g == null || !g.isCalled() || alleleLikelihoodMap == null || alleleLikelihoodMap.isEmpty()) {
            return;
        }
        final Set<Allele> alleles = new HashSet<>(vc.getAlleles());

        // make sure that there's a meaningful relationship between the alleles in the perReadAlleleLikelihoodMap and our VariantContext
        if ( ! alleleLikelihoodMap.getAllelesSet().containsAll(alleles) ) {
            throw new IllegalStateException("VC alleles " + alleles + " not a strict subset of per read allele map alleles " + alleleLikelihoodMap.getAllelesSet());
        }

        final Map<Allele, Integer> alleleCounts = new HashMap<>();
        for ( final Allele allele : vc.getAlleles() ) {
            alleleCounts.put(allele, 0);
        }
        alleleLikelihoodMap.getLikelihoodReadMap().values().stream()
                .map(m -> PerReadAlleleLikelihoodMap.getMostLikelyAllele(m, alleles))
                .filter(a -> a.isInformative())
                .forEach(a -> alleleCounts.compute(a.getMostLikelyAllele(), (allele,prevCount) -> prevCount + 1));

        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference()); //first one in AD is always ref
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++) {
            counts[i + 1] = alleleCounts.get(vc.getAlternateAllele(i));
        }

        gb.AD(counts);
    }

    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.GENOTYPE_ALLELE_DEPTHS); }

    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getFormatLine(getKeyNames().get(0)));
    }
}