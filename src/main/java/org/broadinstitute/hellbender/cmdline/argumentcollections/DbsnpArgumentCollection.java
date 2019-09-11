package org.broadinstitute.hellbender.cmdline.argumentcollections;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureInput;

import java.io.Serializable;

public final class DbsnpArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;
    public static final String DBSNP_LONG_NAME = "dbsnp";
    public static final String DBSNP_SHORT_NAME = "D";

    /**
     * A dbSNP VCF file.
     */
    @Argument(fullName= DBSNP_LONG_NAME, shortName = DBSNP_SHORT_NAME, doc="dbSNP file", optional=true)
    public FeatureInput<VariantContext> dbsnp;

}

