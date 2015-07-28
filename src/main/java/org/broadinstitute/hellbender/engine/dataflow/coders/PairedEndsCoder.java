package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.coders.Coder;
import com.google.cloud.dataflow.sdk.coders.CoderException;
import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.PairedEnds;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Custom coder for the PairedEnds struct-like class.
 */
public class PairedEndsCoder extends CustomCoder<PairedEnds> {
        private static final long serialVersionUID = 1L;

    @Override
    public void encode(PairedEnds value, OutputStream outStream, Context context) throws CoderException, IOException {
        if ( value == null || value.first() == null ) {
            throw new IOException("nothing to encode");
        }

        final boolean isCompletePair = value.second() != null ;
        SerializableCoder.of(Boolean.class).encode(isCompletePair, outStream, context);
        GATKReadCoder gatkReadCoder = new GATKReadCoder();

        gatkReadCoder.encode(value.first(), outStream, context);
        if ( isCompletePair ) {
            gatkReadCoder.encode(value.second(), outStream, context);
        }
    }

    @Override
    public PairedEnds decode(InputStream inStream, Context context) throws CoderException, IOException {
        final boolean isCompletePair = SerializableCoder.of(Boolean.class).decode(inStream, context);
        GATKReadCoder gatkReadCoder = new GATKReadCoder();

        PairedEnds pairedEnds = PairedEnds.of(gatkReadCoder.decode(inStream, context));
        if ( isCompletePair ) {
            pairedEnds.and(gatkReadCoder.decode(inStream, context));
        }
        return pairedEnds;
    }
}
