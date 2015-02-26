package org.broadinstitute.hellbender.cmdline;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureInput;

import java.util.List;

/**
 * An argument collection for use with tools that accept multiple input files containing VariantContext records
 * (eg., VCF/BCF files).
 */
public class StandardVariantInputArgumentCollection implements ArgumentCollectionDefinition {

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more files containing variants", optional = true)
    public List<FeatureInput<VariantContext>> variants;

}
