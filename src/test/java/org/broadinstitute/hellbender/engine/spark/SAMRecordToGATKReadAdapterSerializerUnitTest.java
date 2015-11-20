package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoRegistrator;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.serializer.SerializerInstance;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

public class SAMRecordToGATKReadAdapterSerializerUnitTest {

    public static class TestGATKRegistrator implements KryoRegistrator {
        @SuppressWarnings("unchecked")
        @Override
        public void registerClasses(Kryo kryo) {
            kryo.register(SAMRecordToGATKReadAdapter.class, new SAMRecordToGATKReadAdapterSerializer());
        }
    }

    @Test
    public void testSerializerRoundTripHeaderlessRead() {
        SparkConf conf = new SparkConf().set("spark.kryo.registrator",
                "org.broadinstitute.hellbender.engine.spark.SAMRecordToGATKReadAdapterSerializerUnitTest$TestGATKRegistrator");
        KryoSerializer kryoSerializer = new KryoSerializer(conf);
        SerializerInstance sparkSerializer = kryoSerializer.newInstance();

        // check round trip with no header
        GATKRead read = ArtificialReadUtils.createHeaderlessSamBackedRead("read1", "1", 100, 50);
        check(sparkSerializer, read);
    }

    @Test
    public void testChangingContigsOnHeaderlessGATKRead(){
        final SparkConf conf = new SparkConf().set("spark.kryo.registrator",
                "org.broadinstitute.hellbender.engine.spark.SAMRecordToGATKReadAdapterSerializerUnitTest$TestGATKRegistrator");
        final KryoSerializer kryoSerializer = new KryoSerializer(conf);
        final SerializerInstance sparkSerializer = kryoSerializer.newInstance();
        final GATKRead read = ArtificialReadUtils.createHeaderlessSamBackedRead("read1", "1", 100, 50);
        check(sparkSerializer, read);

        read.setPosition("2", 1);
        check(sparkSerializer, read);
    }

    private void check(SerializerInstance ser, GATKRead read) {
        final ClassTag<GATKRead> tag = ClassTag$.MODULE$.apply(SAMRecordToGATKReadAdapter.class);
        Assert.assertEquals(ser.deserialize(ser.serialize(read, tag), tag), read);
    }
}
