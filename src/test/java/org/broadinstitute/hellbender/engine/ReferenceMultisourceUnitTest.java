package org.broadinstitute.hellbender.engine;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.serializer.SerializerInstance;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

public class ReferenceMultisourceUnitTest extends BaseTest {

    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    @Test
    public void testSerializeRoundTrip2Bit() {
        SparkConf conf = new SparkConf();
        KryoSerializer kryoSerializer = new KryoSerializer(conf);
        SerializerInstance sparkSerializer = kryoSerializer.newInstance();

        PipelineOptions options =null;
        ReferenceMultiSource referenceMultiSource = new ReferenceMultiSource(options, twoBitRefURL, ReferenceWindowFunctions.IDENTITY_FUNCTION);

        final ClassTag<ReferenceMultiSource> tag = ClassTag$.MODULE$.apply(ReferenceMultiSource.class);
        final ReferenceMultiSource actual = sparkSerializer.deserialize(sparkSerializer.serialize(referenceMultiSource, tag), tag);
        Assert.assertEquals(actual.getReferenceSequenceDictionary(null), referenceMultiSource.getReferenceSequenceDictionary(null),
                "\nActual ref: " + actual.getReferenceSequenceDictionary(null) + "\nExpected ref: " + referenceMultiSource.getReferenceSequenceDictionary(null));

    }

}
