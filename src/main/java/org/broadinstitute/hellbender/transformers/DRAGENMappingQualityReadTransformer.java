package org.broadinstitute.hellbender.transformers;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Read transformer intended to replicate DRAGEN behavior for handling mapping qualities.
 *
 * The current behavior maps it naively into the MQ field with no defense against repeat mapping.
 * TODO handle that, make GT deal directly with new annotation
 */
public final class DRAGENMappingQualityReadTransformer implements ReadTransformer {

    private static final long serialVersionUID = 1L;

    public static String EXTENDED_MAPPING_QUALITY_READ_TAG = "XQ";

    private static final int[] mQTableX = new int[]{0, 30, 60, 100, 200, 256};
    private static final int[] mQTableY = new int[]{0, 30, 40, 45, 50, 50};

    public DRAGENMappingQualityReadTransformer( ) { }

    @VisibleForTesting
    public int mapValue(final int val) {
        for (int i = 1; i < mQTableX.length; i++) {
            if (val <= mQTableX[i]) {
                final double xfactor = 1.0*(val - mQTableX[i-1]) / (mQTableX[i] - mQTableX[i-1]) ;
                return  mQTableY[i-1] + (int)(xfactor * (mQTableY[i] - mQTableY[i-1]));
            }
        }
        //This is technically a failure state because the MQ is invalid if somehow it exceeds 256 by the SAMSpec....
        throw new GATKException("Something went wrong trying to map an an invalid XQ tag '"+val+"' val must fall between 0 and 250");
    }

    //TODO figure out how to prevent multiple mappings of a read....
    // Or we just won't worry about it
    @Override
    public GATKRead apply( final GATKRead read ) {
        int oldMQ = read.getMappingQuality();
        int mappedMQ = mapValue(oldMQ);
        read.setMappingQuality(mappedMQ);

//        if (read.hasAttribute(EXTENDED_MAPPING_QUALITY_READ_TAG)) {
//            int extendedMQ = read.getAttributeAsInteger(EXTENDED_MAPPING_QUALITY_READ_TAG);
//            int mappedMQ = mapValue(extendedMQ);
//            read.setMappingQuality(mappedMQ);
//        }
        return read;
    }
}