package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.coders.*;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class ReadContextDataCoder extends CustomCoder<ReadContextData> {
    private static final long serialVersionUID = 1L;
    @Override
    public void encode(ReadContextData value, OutputStream outStream, Context context) throws IOException {
        ReferenceBases overlappingReferenceBases = value.getOverlappingReferenceBases();
        SerializableCoder.of(ReferenceBases.class).encode(overlappingReferenceBases, outStream, context);

        Iterable<Variant> overlappingVariants = value.getOverlappingVariants();
        IterableCoder.of(new VariantCoder()).encode(overlappingVariants, outStream, context);
    }

    @Override
    public ReadContextData decode(InputStream inStream, Context context) throws IOException {
        ReferenceBases referenceBases = SerializableCoder.of(ReferenceBases.class).decode(inStream, context);
        Iterable<Variant> variants = IterableCoder.of(new VariantCoder()).decode(inStream, context);

        return new ReadContextData(referenceBases, variants);
    }
}
