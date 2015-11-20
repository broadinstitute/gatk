package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

public class SAMRecordToGATKReadAdapterSerializer extends Serializer<SAMRecordToGATKReadAdapter> {

    private SAMRecordSparkCodec lazyCodec = new SAMRecordSparkCodec();

    @Override
    public void write(Kryo kryo, Output output, SAMRecordToGATKReadAdapter adapter) {
        SAMRecord record = adapter.getEncapsulatedSamRecord();

        // serialize reference names to avoid having to have a header at read time
        output.writeString(record.getReferenceName());
        output.writeString(record.getMateReferenceName());
        lazyCodec.setOutputStream(output);
        lazyCodec.encode(record);

        // clear indexing bin after encoding to ensure all SAMRecords compare properly
        record.setFlags(record.getFlags());
    }

    @Override
    public SAMRecordToGATKReadAdapter read(Kryo kryo, Input input, Class<SAMRecordToGATKReadAdapter> type) {
        final String referenceName = input.readString();
        final String mateReferenceName = input.readString();
        lazyCodec.setInputStream(input);
        final SAMRecord record = lazyCodec.decode();

        // clear indexing bin after decoding to ensure all SAMRecords compare properly
        record.setFlags(record.getFlags());

        // set reference names (and indexes to null)
        record.setReferenceName(referenceName);
        record.setMateReferenceName(mateReferenceName);

        return SAMRecordToGATKReadAdapter.headerlessReadAdapter(record);
    }
}