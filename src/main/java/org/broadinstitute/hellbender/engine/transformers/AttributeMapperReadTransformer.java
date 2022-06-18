package org.broadinstitute.hellbender.engine.transformers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.ByteMapper;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A read transformer to modify an attribute of a read using a provided lookup table
 *
 * code is based on {@link QualityMapperReadTransformer}
 *
 * see {@link ByteMapper} for lookup table format
 *
 */
public class AttributeMapperReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;

    // parameters
    @Argument(fullName = "attribute-mapping-file")
    public GATKPath attributeMappingFile;

    @Argument(fullName = "attribute-mapping-name")
    public String attributeMappingName;


    // locals
    private ByteMapper mapper;

    @Override
    public GATKRead apply(final GATKRead read) {

        // init?
        if ( mapper == null ) {
            mapper = new ByteMapper(attributeMappingFile);
        }

        // has attribute?
        if ( read.hasAttribute(attributeMappingName) ) {

            // map
            final byte[] src = read.getAttributeAsByteArray(attributeMappingName);
            final byte[] dst = new byte[src.length];

            for (int i = 0; i < src.length; i++) {
                dst[i] = mapper.map(src[i]);
            }

            read.setAttribute(attributeMappingName, dst);
        }

        return read;
    }
}
