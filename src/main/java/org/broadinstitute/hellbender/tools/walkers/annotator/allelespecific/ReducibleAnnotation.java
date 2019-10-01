package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.MergedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.util.stream.Collectors;

/**
 * An interface for annotations that are calculated using raw data across samples, rather than the median (or median of median) of samples values
 *
 */
public interface ReducibleAnnotation<T> extends Annotation {

    String getKeyName();

    String getRawKeyName();

    /**
     * Generate the raw data necessary to calculate the annotation. Raw data is the final endpoint for gVCFs.
     * @param ref the reference context for this annotation
     * @param vc the variant context to annotate
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     */
    default void annotateRaw(final VariantContextBuilder builder,
                             final ReferenceContext ref,
                             final VariantContext vc,
                             final ReadLikelihoods<Allele> likelihoods) {
        final T value = computeRaw(builder, ref, vc, likelihoods);
        setRaw(builder, value);
    }

    default void annotate(final VariantContextBuilder builder,
                          final ReferenceContext ref,
                          final VariantContext vc,
                          final ReadLikelihoods<Allele> likelihoods) {
        final T rawValue = computeRaw(builder, ref, vc, likelihoods);
        finalizeAnnotation(builder, rawValue);
    }

    T computeRaw(final VariantContextBuilder builder, final ReferenceContext ref, final VariantContext vc, final ReadLikelihoods<Allele> likelihoods);

    @SuppressWarnings("unchecked")
    default T getRaw(final VariantContext ctx) {
        return (T) ctx.getAttribute(getRawKeyName());
    }

    default void finalizeAnnotation(final VariantContextBuilder ctx, final T rawValue) {
        ctx.rmAttribute(getRawKeyName());
        ctx.attribute(getKeyName(), toFinalizedAnnotation(rawValue));
    }

    Object toFinalizedAnnotation(final T raw);

    default void setRaw(final VariantContextBuilder builder, final T value) {
         builder.attribute(getRawKeyName(), value);
    }

    T reduceRaws(final VariantContextBuilder builder, final List<T> values);

    default void finalizeAnnotation(final VariantContextBuilder builder, List<VariantContext> vcs) {
        final List<T> raws = vcs.stream()
                .map(this::getRaw)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        if (!raws.isEmpty()) {
            final T reduction = reduceRaws(builder, raws);
            finalizeAnnotation(builder, reduction);
        }
    }

    /**
     * Returns the descriptions used for the VCF INFO meta field corresponding to the annotations raw key.
     * @return A list of VCFInfoHeaderLines corresponding to the raw keys added by this annotaiton
     */
    default List<VCFInfoHeaderLine> getRawDescriptions() {
        final List<VCFInfoHeaderLine> lines = new ArrayList<>(1);
        lines.add(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
        return lines;
    }
}
