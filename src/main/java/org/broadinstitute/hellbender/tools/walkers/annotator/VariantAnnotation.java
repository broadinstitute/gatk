package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.util.Lazy;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Interface of all variant annotations. See also InfoFieldAnnotation and GenotypeFieldAnnotation
 */
public interface VariantAnnotation extends Annotation {

    // is this an INFO or FORMAT annotation
    VCFCompoundHeaderLine.SupportedHeaderLineType annotationType();

    // Return the descriptions used for the VCF INFO or FORMAT meta field.
    default List<VCFCompoundHeaderLine> getDescriptions() {
        return getKeyNames().stream().map(key -> {
            switch (annotationType()) {
                case INFO:
                    return GATKVCFHeaderLines.getInfoLine(key, true);
                case FORMAT:
                    return GATKVCFHeaderLines.getFormatLine(key, true);
                default:
                    throw new IllegalStateException("Unsupported annotation type: " + annotationType());

            }}).collect(Collectors.collectingAndThen(Collectors.toList(), ImmutableList::copyOf));
    }

    /**
     * Return the keys
     */
    List<String> getKeyNames();
}