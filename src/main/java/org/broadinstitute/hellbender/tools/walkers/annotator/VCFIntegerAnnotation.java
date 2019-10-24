package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.vcf.VCFHeaderLineType;

public class VCFIntegerAnnotation implements VCFAnnotation<Integer> {

    public VCFIntegerAnnotation(final VCFAnnotationMeta meta) {
        meta.requiresType(VCFHeaderLineType.Integer);
    }

    @Override
    public VCFAnnotationMeta meta() {
        return null;
    }
}
