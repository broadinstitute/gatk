package org.broadinstitute.hellbender.engine.transformers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.ByteMapper;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

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
    public List<GATKPath> attributeMappingFile;

    @Argument(fullName = "attribute-mapping-name")
    public List<String> attributeMappingName;

    @Argument(fullName = "attribute-shelf-value")
    public byte attributShelfValue;


    // locals
    private Map<String, ByteMapper> mappers;

    @Override
    public GATKRead apply(final GATKRead read) {

        // init?
        if ( mappers == null ) {
            Utils.validate(attributeMappingFile.size() == attributeMappingName.size(), "attribute name and files must have the same length");
            mappers = new LinkedHashMap<>();
            for ( int i = 0 ; i < attributeMappingFile.size() ; i++ ) {
                mappers.put(attributeMappingName.get(i), new ByteMapper(attributeMappingFile.get(i)));
            }
        }

        // has attribute?
        mappers.forEach((attr, byteMapper) -> {
            if ( read.hasAttribute(attr) ) {

                // map
                final boolean isString = isStringAttribute(read, attr);
                final byte[] src = read.getAttributeAsByteArray(attr);
                final byte[] dst = new byte[src.length];

                for (int i = 0; i < src.length; i++) {
                    dst[i] = (byte)(byteMapper.map((byte)(src[i] - attributShelfValue)) + attributShelfValue);
                }

                if ( isString ) {
                    read.setAttribute(attr, new String(dst));
                } else {
                    read.setAttribute(attr, dst);
                }
            }
        });

        return read;
    }

    // at this point, there is no method on GATKRead to determine attribute types
    // TODO: figure out if there is a better way of establishing this, or determine if to add a method to GATKRead
    private boolean isStringAttribute(GATKRead read, String attr) {
        if ( read instanceof SAMRecordToGATKReadAdapter ) {
            return ((SAMRecordToGATKReadAdapter)read).getEncapsulatedSamRecord().getAttribute(attr) instanceof String;
        } else {
            return false;
        }
    }
}
