package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Apply an annotation based on aggregation data from all reads supporting each allele.
 */
public abstract class PerAlleleAnnotation extends GenotypeAnnotation {

    /**
     * Calculate annotations for each allele based on given VariantContext and likelihoods for a given genotype's sample
     * and add the annotations to the GenotypeBuilder.  By default annotations are only calculated for alt alleles but
     * implementations may override the {@code includeRefAllele()} method.  See parent class docs in {@link GenotypeAnnotation}.
     */
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        if ( g == null || likelihoods == null ) {
            return;
        }

        final Map<Allele, List<Integer>> values = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        Utils.stream(likelihoods.bestAllelesBreakingTies(g.getSampleName()))
                .filter(ba -> ba.isInformative() && isUsableRead(ba.read))
                .forEach(ba -> getValueForRead(ba.read, vc).ifPresent(v -> values.get(ba.allele).add(v)));

        final int[] statistics = vc.getAlleles().stream().filter(this::includeAllele).mapToInt(a -> aggregate(values.get(a))).toArray();
        gb.attribute(getVcfKey(), statistics);
    }

    private boolean includeAllele(final Allele allele) {
        return allele.isNonReference() || includeRefAllele();
    }

    // this is false by default but implementations may wish to override
    protected boolean includeRefAllele() { return false; }

    private static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(getVcfKey(), includeRefAllele() ? VCFHeaderLineCount.R : VCFHeaderLineCount.A, VCFHeaderLineType.Integer, getDescription()));
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(getVcfKey()); }

    protected abstract OptionalInt getValueForRead(final GATKRead read, final VariantContext vc);
    protected abstract int aggregate(final List<Integer> values);
    protected abstract String getVcfKey();
    protected abstract String getDescription();
}
