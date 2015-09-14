package org.broadinstitute.hellbender.engine.dataflow.coders;


import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.variant.MinimalVariant;
import org.broadinstitute.hellbender.utils.variant.Variant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Unified coder for Variants.
 *
 * Examines the actual implementation backing the interface, and delegates to the appropriate coder.
 * Encodes the actual type in the byte stream so that the appropriate decoder can be used on the other
 * end.
 */
public final class VariantCoder extends CustomCoder<Variant> {
    private static final long serialVersionUID = 1l;

    @Override
    public void encode( Variant value, OutputStream outStream, Context context ) throws IOException {
        SerializableCoder.of(Class.class).encode(value.getClass(), outStream, context);

        if ( value.getClass() == MinimalVariant.class ) {
            SerializableCoder.of(MinimalVariant.class).encode(((MinimalVariant) value), outStream, context);
        }
        else if (value.getClass() == VariantContextVariantAdapter.class) {
            SerializableCoder.of(VariantContextVariantAdapter.class).encode(((VariantContextVariantAdapter) value), outStream, context);
        } else {
            throw new GATKException("Unknown backing of Variant interface" + value.getClass().getCanonicalName());
        }
    }

    @SuppressWarnings("unchecked")
    @Override
    public Variant decode( InputStream inStream, Context context ) throws IOException {
        final Class<?> clazz = SerializableCoder.of(Class.class).decode(inStream, context);

        if ( clazz == MinimalVariant.class ) {
            return SerializableCoder.of(MinimalVariant.class).decode(inStream, context);
        }
        else if ( clazz == VariantContextVariantAdapter.class) {
            return SerializableCoder.of(VariantContextVariantAdapter.class).decode(inStream, context);
        } else {
            throw new GATKException("Unknown backing of Variant interface" + clazz.getCanonicalName());
        }
    }
}
