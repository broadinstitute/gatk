package org.broadinstitute.hellbender.engine.dataflow.coders;


import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Unified coder for GATKReads.
 *
 * Examines the actual implementation backing the interface, and delegates to the appropriate coder.
 * Encodes the actual type in the byte stream so that the appropriate decoder can be used on the other
 * end.
 *
 * Note that when dealing with a GATKRead PCollection that contains a mix of different backing implementations,
 * you may need to explicitly invoke setCoder(new GATKReadCoder()) on your transform, even if GATKReadCoder
 * has already been registered as the coder for GATKRead, otherwise dataflow may be unable to infer the
 * right coder for your PCollection. This is known to be the case with Create.of(Iterable), in particular.
 */
public final class GATKReadCoder extends CustomCoder<GATKRead> {
    private static final long serialVersionUID = 1l;

    @Override
    public void encode( GATKRead value, OutputStream outStream, Context context ) throws IOException {
        SerializableCoder.of(Class.class).encode(value.getClass(), outStream, context);

        if ( value.getClass() == GoogleGenomicsReadToGATKReadAdapter.class ) {
            GoogleGenomicsReadToGATKReadAdapterCoder.of().encode((GoogleGenomicsReadToGATKReadAdapter)value, outStream, context);
        }
        else if ( value.getClass() == SAMRecordToGATKReadAdapter.class ) {
            SerializableCoder.of(SAMRecordToGATKReadAdapter.class).encode(((SAMRecordToGATKReadAdapter)value), outStream, context);
        }
        else {
            throw new GATKException("Unknown backing implementation for GATKRead: " + value.getClass().getCanonicalName());
        }
    }

    @SuppressWarnings("unchecked")
    @Override
    public GATKRead decode( InputStream inStream, Context context ) throws IOException {
        final Class<?> readClass = SerializableCoder.of(Class.class).decode(inStream, context);

        if ( readClass == GoogleGenomicsReadToGATKReadAdapter.class ) {
            return GoogleGenomicsReadToGATKReadAdapterCoder.of().decode(inStream, context);
        }
        else if ( readClass == SAMRecordToGATKReadAdapter.class ) {
            return SerializableCoder.of(SAMRecordToGATKReadAdapter.class).decode(inStream, context);
        }
        else {
            throw new GATKException("Unknown backing implementation for GATKRead: " + readClass.getCanonicalName());
        }
    }
}
