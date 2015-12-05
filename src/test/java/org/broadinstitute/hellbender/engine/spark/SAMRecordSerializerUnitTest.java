package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoRegistrator;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.serializer.SerializerInstance;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

public class SAMRecordSerializerUnitTest {
    public static class TestGATKRegistrator implements KryoRegistrator {
        @SuppressWarnings("unchecked")
        @Override
        public void registerClasses(Kryo kryo) {
            kryo.register(SAMRecord.class, new SAMRecordSerializer());
        }
    }

    @Test
    public void testSerializerRoundTripHeaderlessSAMRecord() {
        SparkConf conf = new SparkConf().set("spark.kryo.registrator",
                "org.broadinstitute.hellbender.engine.spark.SAMRecordSerializerUnitTest$TestGATKRegistrator");
        KryoSerializer kryoSerializer = new KryoSerializer(conf);
        SerializerInstance sparkSerializer = kryoSerializer.newInstance();

        // check round trip with no header
        final SAMRecord read = ((SAMRecordToGATKReadAdapter)ArtificialReadUtils.createHeaderlessSamBackedRead("read1", "1", 100, 50)).getEncapsulatedSamRecord();
        check(sparkSerializer, read);
    }

    @Test
    public void testChangingContigsOnHeaderlessSAMRecord(){
        final SparkConf conf = new SparkConf().set("spark.kryo.registrator",
                "org.broadinstitute.hellbender.engine.spark.SAMRecordSerializerUnitTest$TestGATKRegistrator");
        final KryoSerializer kryoSerializer = new KryoSerializer(conf);
        final SerializerInstance sparkSerializer = kryoSerializer.newInstance();
        final SAMRecord read = ((SAMRecordToGATKReadAdapter)ArtificialReadUtils.createHeaderlessSamBackedRead("read1", "1", 100, 50)).getEncapsulatedSamRecord();
        check(sparkSerializer, read);

        read.setReferenceName("2");
        read.setAlignmentStart(1);
        check(sparkSerializer, read);
    }

    private void check(SerializerInstance ser, SAMRecord read) {
        final ClassTag<SAMRecord> tag = ClassTag$.MODULE$.apply(SAMRecord.class);
        final SAMRecord actual = ser.deserialize(ser.serialize(read, tag), tag);
        Assert.assertEquals(actual, read, "\nActual read: " + actual.getSAMString() + "\nExpected read: " + read.getSAMString());
    }
}
