package org.broadinstitute.hellbender.cmdline.argumentcollections;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureInput;

import java.io.Serializable;
import java.util.List;

/**
 * An argument collection for use with tools that accept one or more input files containing VariantContext records
 * (eg., VCF files), and require at least one such input.
 */
public final class RequiredVariantInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more files containing variants", optional = false)
    public List<FeatureInput<VariantContext>> variantFiles;

}
