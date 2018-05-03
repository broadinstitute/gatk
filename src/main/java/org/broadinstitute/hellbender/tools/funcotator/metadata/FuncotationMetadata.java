package org.broadinstitute.hellbender.tools.funcotator.metadata;

import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.List;

/**
 * Represents metadata information for fields in in a Funcotation.
 *
 * Even though types are from the VCF Spec, this can be used for other formats besides VCF.
 *
 */
public interface FuncotationMetadata {

    /**
     * @param fieldName field to search.  Never {@code null}
     * @return Returns the metadata for the given field name as a VCF Info Header Line.  {@code null} is returned if
     * the given field name is not in this metadata.
     */
    VCFInfoHeaderLine retrieveHeaderInfo(final String fieldName);

    /**
     * @return All of the fields in this metadata object.  Never {@code null}.  Empty list if this metadata instance has
     *  no fields.
     */
    List<VCFInfoHeaderLine> retrieveAllHeaderInfo();
}
