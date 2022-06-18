package org.broadinstitute.hellbender.engine.transformers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.ByteMapper;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * A read transformer to modify the quality of a read using a provided lookup table
 *
 * see {@link ByteMapper} for lookup table format
 */
public class QualityMapperReadTransformer implements ReadTransformer {
    private static final long serialVersionUID = 1L;

    // parameters
    @Argument(fullName = "quality-mapping-file")
    public GATKPath qualityMappingFile;

    // locals
    private ByteMapper mapper;

    @Override
    public GATKRead apply(final GATKRead read) {

        // init?
        if ( mapper == null ) {
            mapper = new ByteMapper(qualityMappingFile);
        }

        // map
        final byte[] qual = read.getBaseQualitiesNoCopy();
        final byte[] newQual = new byte[qual.length];

        for ( int i = 0 ; i < qual.length ; i++ ) {
            newQual[i] = mapper.map(qual[i]);
        }

        read.setBaseQualities(newQual);

        return read;
    }
}
