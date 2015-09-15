package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

public class SAMRecordToGATKReadAdapterSerializer extends Serializer<SAMRecordToGATKReadAdapter> {

    private BAMRecordCodec lazyCodec = new BAMRecordCodec(null, new LazyBAMRecordFactory());

    @Override
    public void write(Kryo kryo, Output output, SAMRecordToGATKReadAdapter adapter) {
        SAMRecord record = adapter.getSamRecord();
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

        final int referenceIndex = record.getReferenceIndex();
        final int mateReferenceIndex = record.getMateReferenceIndex();
        record.setReferenceName(referenceName);
        record.setMateReferenceName(mateReferenceName);

        // set reference indexes again since setting reference names without a header will unset indexes
        record.setReferenceIndex(referenceIndex);
        record.setMateReferenceIndex(mateReferenceIndex);

        return (SAMRecordToGATKReadAdapter) SAMRecordToGATKReadAdapter.sparkReadAdapter(record);
    }


    static class LazyBAMRecordFactory implements SAMRecordFactory {
        @Override public SAMRecord createSAMRecord(SAMFileHeader hdr) {
            throw new UnsupportedOperationException(
                    "LazyBAMRecordFactory can only create BAM records");
        }

        @Override public BAMRecord createBAMRecord(
                SAMFileHeader hdr,
                int referenceSequenceIndex, int alignmentStart,
                short readNameLength, short mappingQuality,
                int indexingBin, int cigarLen, int flags, int readLen,
                int mateReferenceSequenceIndex, int mateAlignmentStart,
                int insertSize, byte[] variableLengthBlock)
        {
            return new LazyBAMRecord(
                    hdr, referenceSequenceIndex, alignmentStart, readNameLength,
                    mappingQuality, indexingBin, cigarLen, flags, readLen,
                    mateReferenceSequenceIndex, mateAlignmentStart, insertSize,
                    variableLengthBlock);
        }
    }

    static class LazyBAMRecord extends BAMRecord {
        private static final long serialVersionUID = 1L;

        public LazyBAMRecord(
                SAMFileHeader hdr, int referenceID, int coordinate, short readNameLength,
                short mappingQuality, int indexingBin, int cigarLen, int flags,
                int readLen, int mateReferenceID, int mateCoordinate, int insertSize,
                byte[] restOfData)
        {
            super(
                    hdr, referenceID, coordinate, readNameLength, mappingQuality,
                    indexingBin, cigarLen, flags, readLen, mateReferenceID,
                    mateCoordinate, insertSize, restOfData);
        }

        @Override
        public void setReferenceIndex(final int referenceIndex) {
            mReferenceIndex = referenceIndex;
        }

        @Override
        public void setMateReferenceIndex(final int referenceIndex) {
            mMateReferenceIndex = referenceIndex;
        }
    }
}