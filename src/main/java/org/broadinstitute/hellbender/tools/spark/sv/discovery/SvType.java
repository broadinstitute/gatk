package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.variant.variantcontext.Allele;

import java.util.Collections;
import java.util.Map;


/**
 * Various types of structural variations
 */
public abstract class SvType {

    protected final String variantId;
    protected final Allele altAllele;
    protected final int svLen;
    protected final Map<String, String> extraAttributes;

    protected static final Map<String, String> noExtraAttributes = Collections.emptyMap();

    protected SvType(final String id, final Allele altAllele, final int len, final Map<String, String> typeSpecificExtraAttributes) {
        variantId = id;
        this.altAllele = altAllele;
        svLen = len;
        extraAttributes = typeSpecificExtraAttributes;
    }

    public final String getInternalVariantId() {
        return variantId;
    }
    public final Allele getAltAllele() {
        return altAllele;
    }
    public final int getSVLength() {
        return svLen;
    }
    public final Map<String, String> getTypeSpecificAttributes() {
        return extraAttributes;
    }
    public abstract boolean isBreakEndOnly();
}
