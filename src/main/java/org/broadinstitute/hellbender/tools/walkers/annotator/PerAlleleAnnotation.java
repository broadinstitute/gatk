package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Apply an annotation based on aggregation data from all reads supporting each allele.
 */
public abstract class PerAlleleAnnotation extends InfoFieldAnnotation{

    /**
     * Calculate annotations for each allele based on given VariantContext and likelihoods for a given genotype's sample
     * and add the annotations to the GenotypeBuilder.  By default annotations are only calculated for alt alleles but
     * implementations may override the {@code includeRefAllele()} method.  See parent class docs in {@link GenotypeAnnotation}.
     */
    public Map<String, Object> annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        if ( likelihoods == null ) {
            return Collections.emptyMap();
        }

        final Map<Allele, List<Integer>> values = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        Utils.stream(likelihoods.bestAllelesBreakingTies())
                .filter(ba -> ba.isInformative() && isUsableRead(ba.evidence))
                .forEach(ba -> getValueForRead(ba.evidence, vc).ifPresent(v -> values.get(ba.allele).add(v)));

        final int[] statistics = vc.getAlleles().stream().filter(this::includeAllele).mapToInt(a -> aggregate(values.get(a))).toArray();
        return ImmutableMap.of(getVcfKey(), statistics);
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
    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine(getVcfKey(), includeRefAllele() ? VCFHeaderLineCount.R : VCFHeaderLineCount.A, VCFHeaderLineType.Integer, getDescription()));
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(getVcfKey()); }

    protected abstract OptionalInt getValueForRead(final GATKRead read, final VariantContext vc);
    protected abstract int aggregate(final List<Integer> values);
    protected abstract String getVcfKey();
    protected abstract String getDescription();
}
