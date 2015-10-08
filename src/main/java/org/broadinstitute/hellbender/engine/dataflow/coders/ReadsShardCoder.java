package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.coders.CustomCoder;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.UUID;

/**
 * A coder for ContextShard that is more efficient than Java serialization.
 * The main trick is to use the BAMRecordCodec for the reads.
 * It assumes the reads are all SAMRecordToGATKReadAdapter.
 */
public final class ReadsShardCoder extends CustomCoder<ReadsShard> {
    private static final long serialVersionUID = 1l;
    private transient BAMRecordCodec lazyCodec = null;
    private final byte[] buffer = new byte[32];

    private synchronized void init() {
        // working around BAMRecordCodec not being serializable
        if (null==lazyCodec) {
            lazyCodec = new BAMRecordCodec(null, new org.broadinstitute.hellbender.engine.spark.SAMRecordToGATKReadAdapterSerializer.LazyBAMRecordFactory());
        }
    }

    @Override
    public void encode(ReadsShard value, OutputStream outStream, Context context) throws IOException {
        init();
        writeSimpleInterval(outStream, value.interval);
        lazyCodec.setOutputStream(outStream);
        writeInt(outStream, value.reads.size());
        String refName = null;
        String mateRefName = null;
        for (GATKRead g : value.reads) {
            writeUUID(outStream, g.getUUID());
            SAMRecord r = ((SAMRecordToGATKReadAdapter)g).convertToSAMRecord(null);
            lazyCodec.encode(r);
            String s = r.getReferenceName();
            boolean writeRefName = (!s.equals(refName));
            boolean mateHasRefName = r.getMateReferenceIndex()!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
            boolean writeMateRefName = mateHasRefName;
            String t = null;
            if (mateHasRefName) {
                t = r.getMateReferenceName();
                if (t.equals(mateRefName)) writeMateRefName=false;
            }
            // not writing the reference name when it's the same saves us about 2% space and a bit more time.

            // 000 = same refName, no mateRefName
            // 001 = same refName, same mateRefName
            // 011 = same refName, new mateRefName
            // 1.. = new refName, ...
            byte indicator = 0;
            indicator += (writeRefName?4:0);
            indicator += (writeMateRefName?2:0);
            indicator += (mateHasRefName?1:0);

            outStream.write(indicator);
            if (writeRefName) {
                refName = s;
                writeString(outStream, refName);
            }
            if (writeMateRefName) {
                mateRefName = t;
                writeString(outStream, mateRefName);
            }
            // clear indexing bin after encoding to ensure all SAMRecords compare properly
            r.setFlags(r.getFlags());
        }
    }

    @Override
    public ReadsShard decode(InputStream inStream, Context context) throws IOException {
        init();
        SimpleInterval interval = readSimpleInterval(inStream);
        String refName = null;
        String mateRefName = null;
        lazyCodec.setInputStream(inStream);
        int size = readInt(inStream);
        ArrayList<GATKRead> reads = new ArrayList<>(size);
        for (int i=0; i<size; i++) {
            UUID uuid = readUUID(inStream);
            SAMRecord record = lazyCodec.decode();
            final int referenceIndex = record.getReferenceIndex();
            final int mateReferenceIndex = record.getMateReferenceIndex();

            int indicator = inStream.read();
            if (indicator<0) throw new EOFException("EOF partway through shard record");

            boolean mateHasRefName = (indicator&1)!=0;
            boolean readMateRefName = (indicator&2)!=0;
            boolean readRefName = (indicator&4)!=0;

            if (readRefName) {
                refName = readString(inStream);
            }
            record.setReferenceName(refName);
            if (mateHasRefName) {
                if (readMateRefName) {
                    mateRefName = readString(inStream);
                }
                record.setMateReferenceName(mateRefName);
            }

            // clear indexing bin after decoding to ensure all SAMRecords compare properly
            record.setFlags(record.getFlags());

            // set reference indexes again since setting reference names without a header will unset indexes
            record.setReferenceIndex(referenceIndex);
            record.setMateReferenceIndex(mateReferenceIndex);

            reads.add(new SAMRecordToGATKReadAdapter(record, uuid));

        }
        return new ReadsShard(interval, reads);
    }

    private void writeSimpleInterval(OutputStream o, SimpleInterval i) throws IOException {
        writeString(o, i.getContig());
        writeInt(o, i.getStart());
        writeInt(o, i.getEnd());
    }

    private SimpleInterval readSimpleInterval(InputStream i) throws IOException {
        String contig = readString(i);
        int start = readInt(i);
        int end = readInt(i);
        return new SimpleInterval(contig, start, end);
    }

    private void writeUUID(OutputStream o, UUID id) throws IOException {
        ByteBuffer.wrap(buffer).putLong(id.getLeastSignificantBits()).putLong(id.getMostSignificantBits());
        o.write(buffer, 0, 16);
    }

    private UUID readUUID(InputStream i) throws IOException {
        i.read(buffer, 0, 16);
        ByteBuffer buf = ByteBuffer.wrap(buffer);
        long lsb = buf.getLong();
        long msb = buf.getLong();
        return new UUID(msb, lsb);
    }

    private void writeInt(OutputStream o, int i) throws IOException {
        ByteBuffer.wrap(buffer).putInt(i).array();
        o.write(buffer, 0, 4);
    }

    private int readInt(InputStream i) throws IOException {
        i.read(buffer, 0, 4);
        return ByteBuffer.wrap(buffer).getInt();
    }

    private void writeString(OutputStream o, String s) throws IOException {
        byte[] b = s.getBytes(StandardCharsets.UTF_8);
        writeInt(o, b.length);
        o.write(b);
    }

    private String readString(InputStream i) throws IOException {
        int len = readInt(i);
        byte[] b = new byte[len];
        i.read(b);
        return new String(b, "UTF-8");
    }
}