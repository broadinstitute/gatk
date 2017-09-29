package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.math3.util.DoubleArray;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Apply an annotation based on aggregation data from all reads supporting each allele.
 *
 * Created by David Benjamin on 4/13/17.
 */
public abstract class PerAlleleAnnotation extends GenotypeAnnotation {

    /**
     * Calculate annotations for eah allele based on given VariantContext and likelihoods for a given genotype's sample
     * and add the annotations to the GenotypeBuilder.  See parent class docs in {@link GenotypeAnnotation}.
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

        Utils.stream(likelihoods.bestAlleles(g.getSampleName()))
                .filter(ba -> ba.isInformative() && isUsableRead(ba.read))
                .forEach(ba -> getValueForRead(ba.read, vc).ifPresent(v -> values.get(ba.allele).add(v)));

        final int[] statistics = vc.getAlleles().stream().mapToInt(a -> aggregate(values.get(a))).toArray();
        gb.attribute(getVcfKey(), statistics);
    }

    static boolean isUsableRead(final GATKRead read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(getVcfKey(), VCFHeaderLineCount.A, VCFHeaderLineType.Float, getDescription()));
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(getVcfKey()); }

    protected abstract OptionalInt getValueForRead(final GATKRead read, final VariantContext vc);
    protected abstract int aggregate(final List<Integer> values);
    protected abstract String getVcfKey();
    protected abstract String getDescription();
}
