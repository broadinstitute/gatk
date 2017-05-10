package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.Map;

/**
 * A template name and an intervalId.
 */
@DefaultSerializer(QNameAndInterval.Serializer.class)
public final class QNameAndInterval implements Map.Entry<String, Integer> {
    private final byte[] qName;
    private final int hashVal;
    private final int intervalId;

    public QNameAndInterval( final String qName, final int intervalId ) {
        this.qName = qName.getBytes();
        this.hashVal = qName.hashCode() ^ (47 * intervalId);
        this.intervalId = intervalId;
    }

    private QNameAndInterval( final Kryo kryo, final Input input ) {
        final int nameLen = input.readInt();
        qName = input.readBytes(nameLen);
        hashVal = input.readInt();
        intervalId = input.readInt();
    }

    private void serialize( final Kryo kryo, final Output output ) {
        output.writeInt(qName.length);
        output.writeBytes(qName);
        output.writeInt(hashVal);
        output.writeInt(intervalId);
    }

    @Override
    public String getKey() {
        return getQName();
    }

    @Override
    public Integer getValue() {
        return getIntervalId();
    }

    @Override
    public Integer setValue( final Integer value ) {
        throw new UnsupportedOperationException("Can't set QNameAndInterval.intervalId");
    }

    public String getQName() {
        return new String(qName);
    }

    public int getIntervalId() {
        return intervalId;
    }

    @Override
    public int hashCode() {
        return hashVal;
    }

    @Override
    public boolean equals( final Object obj ) {
        return obj instanceof QNameAndInterval && equals((QNameAndInterval) obj);
    }

    public boolean equals( final QNameAndInterval that ) {
        return this.intervalId == that.intervalId && Arrays.equals(this.qName, that.qName);
    }

    public String toString() {
        return new String(qName) + " " + intervalId;
    }

    /**
     * write template names and interval IDs to a file.
     */
    public static void writeQNames( final String qNameFile,
                                    final Iterable<QNameAndInterval> qNames ) {
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(qNameFile))) ) {
            for ( final QNameAndInterval qnameAndInterval : qNames ) {
                writer.write(qnameAndInterval.toString() + "\n");
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write qname intervals file " + qNameFile, ioe);
        }
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<QNameAndInterval> {
        @Override
        public void write( final Kryo kryo, final Output output, final QNameAndInterval qNameAndInterval ) {
            qNameAndInterval.serialize(kryo, output);
        }

        @Override
        public QNameAndInterval read( final Kryo kryo, final Input input, final Class<QNameAndInterval> klass ) {
            return new QNameAndInterval(kryo, input);
        }
    }
}
