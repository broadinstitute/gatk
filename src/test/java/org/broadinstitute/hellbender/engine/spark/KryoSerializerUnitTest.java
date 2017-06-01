package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoRegistrator;
import org.broadinstitute.hellbender.utils.collections.IndexedSet;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.test.SparkTestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class KryoSerializerUnitTest {

    @Test
    public void testSerializerRoundTrip() {
        SparkConf conf = new SparkConf();
        IndexedSet<String> set = new IndexedSet<>();
        set.add("a");

        IndexedSet<String> roundTrippedSet = SparkTestUtils.roundTripInKryo(set, IndexedSet.class, conf);
        Assert.assertEquals(roundTrippedSet, set);
    }

}
