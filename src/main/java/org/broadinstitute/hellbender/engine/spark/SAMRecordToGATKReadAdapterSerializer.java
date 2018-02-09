package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;

/**
 * Efficient serializer for SAMRecordToGATKReadAdapters that uses SAMRecordSparkCodec for encoding/decoding.
 * Assumes that the underlying SAMRecords are headerless (and clears their header if they're not).
 */
public final class SAMRecordToGATKReadAdapterSerializer extends Serializer<SAMRecordToGATKReadAdapter> {

    private SAMRecordSparkCodec lazyCodec = new SAMRecordSparkCodec();
    private Map<String, String> stringInternedMap = new HashMap<>();

    final Field mReferenceName;
    final Field mReferenceIndex;
    final Field mMateReferenceName;
    final Field mMateReferenceIndex;

    public SAMRecordToGATKReadAdapterSerializer(){
        super();
        try {
            mReferenceName = SAMRecord.class.getDeclaredField("mReferenceName");
            mReferenceName.setAccessible(true);
            mReferenceIndex = SAMRecord.class.getDeclaredField("mReferenceIndex");
            mReferenceIndex.setAccessible(true);

            mMateReferenceName = SAMRecord.class.getDeclaredField("mMateReferenceName");
            mMateReferenceName.setAccessible(true);
            mMateReferenceIndex = SAMRecord.class.getDeclaredField("mMateReferenceIndex");
            mMateReferenceIndex.setAccessible(true);
        } catch (NoSuchFieldException e) {
            throw new RuntimeException(e);
        }

    }

    @Override
    public void write(Kryo kryo, Output output, SAMRecordToGATKReadAdapter adapter) {
        SAMRecord record = adapter.getEncapsulatedSamRecord();
        // The underlying read is likely to already be headerless, but as a defensive
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
    public SAMRecordToGATKReadAdapter read(Kryo kryo, Input input, Class<SAMRecordToGATKReadAdapter> type) {
        final String referenceName = input.readString();
        final String mateReferenceName = input.readString();
        lazyCodec.setInputStream(input);
        final SAMRecord record = lazyCodec.decode();

        // clear indexing bin after decoding to ensure all SAMRecords compare properly
        record.setFlags(record.getFlags());
        
        try {
            if (SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(referenceName)) {
                mReferenceName.set(record, SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
                mReferenceIndex.set(record, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            } else {
                mReferenceName.set(record, stringInternedMap.computeIfAbsent(referenceName, String::intern));
                mReferenceIndex.set(record, null);
            }
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }

        try {
            if (SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(mateReferenceName)) {
                mMateReferenceName.set(record, SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
                mMateReferenceIndex.set(record, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            } else {
                mMateReferenceName.set(record, stringInternedMap.computeIfAbsent(mateReferenceName, String::intern));
                mMateReferenceIndex.set(record, null);
            }
        } catch (IllegalAccessException e) {
            throw new RuntimeException(e);
        }

        // headerlessReadAdapter() calls setHeaderStrict(null), which will set reference indices to null if the above
        // setReferenceName()/setMateReferenceName() calls failed to do so (eg., in the case of "*" as the
        // reference name).
        return SAMRecordToGATKReadAdapter.headerlessReadAdapter(record);
    }
}