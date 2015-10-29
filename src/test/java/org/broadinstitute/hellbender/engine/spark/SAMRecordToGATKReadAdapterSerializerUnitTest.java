package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
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
    public void test() {
        SparkConf conf = new SparkConf().set("spark.kryo.registrator",
                "org.broadinstitute.hellbender.engine.spark.SAMRecordToGATKReadAdapterSerializerUnitTest$TestGATKRegistrator");
        KryoSerializer kryoSerializer = new KryoSerializer(conf);
        SerializerInstance sparkSerializer = kryoSerializer.newInstance();

        // check round trip with header set
        GATKRead read = ArtificialReadUtils.createSamBackedRead("read1", "1", 100, 50);
        check(sparkSerializer, read);

        // check round trip with no header
        ((SAMRecordToGATKReadAdapter) read).getSamRecord().setHeader(null);
        check(sparkSerializer, read);
    }

    private void check(SerializerInstance ser, GATKRead read) {
        final ClassTag<GATKRead> tag = ClassTag$.MODULE$.apply(SAMRecordToGATKReadAdapter.class);
        Assert.assertEquals(ser.deserialize(ser.serialize(read, tag), tag), read);
    }
}
