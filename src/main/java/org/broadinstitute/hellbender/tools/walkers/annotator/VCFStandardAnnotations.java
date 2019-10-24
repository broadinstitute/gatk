package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

public final class VCFStandardAnnotations {

    public static final VCFIntegerAnnotation VAR_DP = VCFAnnotations.builder().target(VCFAnnotationTarget.VARIANT)
                                                                    .type(VCFHeaderLineType.Integer)
                                                                    .number(1)
                                                                    .description("Depth of a site in reads or fragments")
                                                                    .make(VCFIntegerAnnotation::new);

}
