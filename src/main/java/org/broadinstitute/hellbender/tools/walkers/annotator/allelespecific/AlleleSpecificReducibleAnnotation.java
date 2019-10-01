package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public interface AlleleSpecificReducibleAnnotation<T> extends ReducibleAnnotation<Map<Allele, T>> {

    T reduce(final Allele allele, List<T> alleleValues);

    @Override
    default Map<Allele, T> reduce(final List<Map<Allele, T>> values) {
        if (values.isEmpty()) {
            return Collections.emptyMap();
        } else {
            final Map<Allele, List<T>> listed = values.stream()
                    .flatMap(value -> value.entrySet().stream())
                    .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(Map.Entry::getValue, Collectors.toList())));
            final Map<Allele, T> result = new HashMap<>(listed.size());
            for (final Map.Entry<Allele, List<T>> entry : listed.entrySet()) {
                result.put(entry.getKey(), reduce(entry.getKey(), entry.getValue()));
            }
            return result;
        }
    }

    T computeRaw(final int alleleIndex, final Allele allele, final ReferenceContext ref, final VariantContext ctx, final ReadLikelihoods<Allele> likelihoods);

    @Override
    default Map<Allele, T> computeRaw(final ReferenceContext ref, final VariantContext vc, final ReadLikelihoods<Allele> likelihoods) {
        final List<Allele> alleles = vc.getAlleles();
        final Map<Allele, T> result = new HashMap<>();
        for (int i = 0; i < alleles.size(); i++) {
            final Allele allele = alleles.get(0);
            final T value = computeRaw(i, allele, ref, vc, likelihoods);
            result.put(allele, value);
        }
        return result;
    }
}
