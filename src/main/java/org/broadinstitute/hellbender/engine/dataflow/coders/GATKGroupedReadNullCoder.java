package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.values.KV;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;

/**
 * This coder does nothing; we can use it in fused stages to reduce sampling overheads.
 */
public final class GATKGroupedReadNullCoder extends CustomCoder<KV<String, Iterable<GATKRead>>> {
    private static final long serialVersionUID = 1l;
    private boolean allowDecode = false;

    public GATKGroupedReadNullCoder() {}

    public GATKGroupedReadNullCoder(boolean allowDecode) {
        this.allowDecode = allowDecode;
    }


    @Override
    public void encode( KV<String, Iterable<GATKRead>> value, OutputStream outStream, Context context ) throws IOException {
    }

    @SuppressWarnings("unchecked")
    @Override
    public KV<String, Iterable<GATKRead>> decode( InputStream inStream, Context context ) throws IOException {
        if (!allowDecode) {
            throw new RuntimeException("cannot decode the null stream.");
        }
        return KV.of("null-coder", new ArrayList<GATKRead>());
    }
}
