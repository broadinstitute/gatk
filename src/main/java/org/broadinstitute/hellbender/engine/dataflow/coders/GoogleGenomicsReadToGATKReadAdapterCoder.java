package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.coders.CustomCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.genomics.dataflow.coders.GenericJsonCoder;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.UUID;

public class GoogleGenomicsReadToGATKReadAdapterCoder extends CustomCoder<GoogleGenomicsReadToGATKReadAdapter> {
    private static final long serialVersionUID = 1l;

    public static GoogleGenomicsReadToGATKReadAdapterCoder of() {
        return new GoogleGenomicsReadToGATKReadAdapterCoder();
    }

    @Override
    public String getEncodingId() {
        return "GoogleGenomicsReadToGATKReadAdapterCoder";
    }

    @Override
    public void encode(GoogleGenomicsReadToGATKReadAdapter value, OutputStream outStream, Context context) throws IOException {
        KvCoder.of(UUIDCoder.CODER, GenericJsonCoder.of(Read.class)).encode(KV.of(value.getUUID(), value.convertToGoogleGenomicsRead()), outStream, context);
    }

    @Override
    public GoogleGenomicsReadToGATKReadAdapter decode(InputStream inStream, Context context) throws IOException {
        final KV<UUID, Read> decode = KvCoder.of(UUIDCoder.CODER, GenericJsonCoder.of(Read.class)).decode(inStream, context);
        final UUID uuid = decode.getKey();
        final Read read = decode.getValue();
        return new GoogleGenomicsReadToGATKReadAdapter(read, uuid);
    }
}
