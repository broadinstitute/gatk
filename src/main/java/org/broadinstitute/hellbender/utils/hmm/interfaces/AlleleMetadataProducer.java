package org.broadinstitute.hellbender.utils.hmm.interfaces;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeader;

import javax.annotation.Nonnull;

/**
 * An interface for hidden states which are construable as {@link Allele}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public interface AlleleMetadataProducer {

    /**
     * Get the corresponding {@link Allele}
     * @return an instance of {@link Allele}
     */
    Allele toAllele();

    /**
     * Append relevant header lines to a {@link VCFHeader} instance
     * @param header an instance of {@link VCFHeader}
     */
    void addHeaderLineTo(@Nonnull final VCFHeader header);
}
