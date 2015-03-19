package org.broadinstitute.hellbender.cmdline.argumentcollections;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.engine.FeatureInput;

public final class DbsnpArgumentCollection implements ArgumentCollectionDefinition{

    /**
     * A dbSNP VCF file.
     */
    @Argument(fullName="dbsnp", shortName = "D", doc="dbSNP file", optional=true)
    public FeatureInput<VariantContext> dbsnp;

}

