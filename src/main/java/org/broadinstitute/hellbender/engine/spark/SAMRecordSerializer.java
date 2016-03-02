package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.*;

/**
 * Efficient serializer for SAMRecords that uses SAMRecordSparkCodec for encoding/decoding.
 * Assumes that the SAMRecords are headerless (and clears their header if they're not).
 */
public final class SAMRecordSerializer extends Serializer<SAMRecord> {
    private SAMRecordSparkCodec lazyCodec = new SAMRecordSparkCodec();

    @Override
    public void write(Kryo kryo, Output output, SAMRecord record) {
        // The read is likely to already be headerless, but as a defensive
        // measure in case it's not, set the header to null explicitly.
        record.setHeaderStrict(null);

        // serialize reference names to avoid having to have a header at read time
        output.writeString(record.getReferenceName());
        output.writeString(record.getMateReferenceName());
        lazyCodec.setOutputStream(output);
        lazyCodec.encode(record);

        // clear indexing bin after encoding to ensure all SAMRecords compare properly
        record.setFlags(record.getFlags());
    }

    @Override
    public SAMRecord read(Kryo kryo, Input input, Class<SAMRecord> type) {
        final String referenceName = input.readString();
        final String mateReferenceName = input.readString();
        lazyCodec.setInputStream(input);
        final SAMRecord record = lazyCodec.decode();

        // clear indexing bin after decoding to ensure all SAMRecords compare properly
        record.setFlags(record.getFlags());

        // set reference names (and indices to null)
        record.setReferenceName(referenceName);
        record.setMateReferenceName(mateReferenceName);
        // Explicitly clear the reference indices by calling setHeaderStrict(null). Although setReferenceName()
        // and setMateReferenceName() above will usually null out the reference indices for us (since our
        // read is headerless) they won't do so if either name is "*"
        record.setHeaderStrict(null);

        return record;
    }
}
