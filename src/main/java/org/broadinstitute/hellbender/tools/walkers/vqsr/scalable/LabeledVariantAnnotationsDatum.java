package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.List;
import java.util.Set;

final class LabeledVariantAnnotationsDatum implements Locatable {
    public final SimpleInterval loc;
    public final Allele refAllele;
    public final Allele altAllele;
    public final ImmutableSet<String> labels;
    public final ImmutableList<Double> annotations;
    public double score;    // TODO decouple

    public LabeledVariantAnnotationsDatum(final VariantContext vc,
                                          final Allele refAllele,
                                          final Allele altAllele,
                                          final Set<String> labels,
                                          final List<String> sortedAnnotationKeys,
                                          final boolean useASAnnotations ) {
        this.loc = new SimpleInterval(vc);
        this.refAllele = refAllele;
        this.altAllele = altAllele;
        this.labels = ImmutableSet.copyOf(labels);
        this.annotations = sortedAnnotationKeys.stream()
                .map(k -> decodeAnnotation(k, vc, altAllele, useASAnnotations))
                .collect(ImmutableList.toImmutableList());
    }

    @Override
    public String getContig() {
        return loc.getContig();
    }

    @Override
    public int getStart() {
        return loc.getStart();
    }

    @Override
    public int getEnd() {
        return loc.getEnd();
    }

    private static double decodeAnnotation(final String annotationKey,
                                           final VariantContext vc,
                                           final Allele altAllele,
                                           final boolean useASannotations) {
        double value;

        try {
            //if we're in allele-specific mode and an allele-specific annotation has been requested, parse the appropriate value from the list
            if (useASannotations && annotationKey.startsWith(GATKVCFConstants.ALLELE_SPECIFIC_PREFIX)) {
                final List<Object> valueList = vc.getAttributeAsList(annotationKey);
                //FIXME: we need to look at the ref allele here too
                if (vc.hasAllele(altAllele)) {
                    final int altIndex = vc.getAlleleIndex(altAllele) - 1; //- 1 is to convert the index from all alleles (including reference) to just alternate alleles
                    value = Double.parseDouble((String) valueList.get(altIndex));
                } else {
                    //if somehow our alleles got mixed up
                    throw new IllegalStateException("ExtractAnnotationsVariantDatum allele " + altAllele + " is not contained in the input VariantContext.");
                }
            } else {
                value = vc.getAttributeAsDouble(annotationKey, Double.NaN);
            }
            if (Double.isInfinite(value)) {
                value = Double.NaN;
            }
        } catch (final NumberFormatException e) {
            value = Double.NaN;
        }
        return value;
    }
}