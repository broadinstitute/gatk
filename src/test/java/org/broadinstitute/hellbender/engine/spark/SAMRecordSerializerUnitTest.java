package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoRegistrator;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.testutils.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

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

        // check round trip with no header
        final SAMRecord read = ((SAMRecordToGATKReadAdapter)ArtificialReadUtils.createHeaderlessSamBackedRead("read1", "1", 100, 50)).getEncapsulatedSamRecord();

        final SAMRecord roundTrippedRead = SparkTestUtils.roundTripInKryo(read, SAMRecord.class, conf);
        Assert.assertEquals(roundTrippedRead, read, "\nActual read: " + roundTrippedRead.getSAMString() + "\nExpected read: " + read.getSAMString());
    }

    @Test
    public void testChangingContigsOnHeaderlessSAMRecord(){
        final SparkConf conf = new SparkConf().set("spark.kryo.registrator",
                "org.broadinstitute.hellbender.engine.spark.SAMRecordSerializerUnitTest$TestGATKRegistrator");
        final SAMRecord read = ((SAMRecordToGATKReadAdapter)ArtificialReadUtils.createHeaderlessSamBackedRead("read1", "1", 100, 50)).getEncapsulatedSamRecord();

        final SAMRecord roundTrippedRead = SparkTestUtils.roundTripInKryo(read, SAMRecord.class, conf);
        Assert.assertEquals(roundTrippedRead, read, "\nActual read: " + roundTrippedRead.getSAMString() + "\nExpected read: " + read.getSAMString());

        read.setReferenceName("2");
        read.setAlignmentStart(1);

        final SAMRecord roundTrippedRead2 = SparkTestUtils.roundTripInKryo(read, SAMRecord.class, conf);
        Assert.assertEquals(roundTrippedRead2, read, "\nActual read: " + roundTrippedRead2.getSAMString() + "\nExpected read: " + read.getSAMString());
    }
}
