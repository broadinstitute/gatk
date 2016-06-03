package org.broadinstitute.hellbender.utils;


import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.serializer.SerializerInstance;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

public class SerializationTestUtils {

    /**
     *  Test that an instance can be serialized in kryo.
     *
     * @param input instance of inputClazz.  Never {@code null}
     * @param inputClazz class to cast input
     * @param conf Spark configuration to test
     * @param <T> class to attempt.  Same or subclass of inputClazz
     * @return serialized and deserialized instance of input.  Throws exception if serialization round trip fails.
     */
    public static <T> T roundTripInKryo(T input, Class<?> inputClazz, final SparkConf conf) {
        Utils.nonNull(input);
        KryoSerializer kryoSerializer = new KryoSerializer(conf);
        SerializerInstance sparkSerializer = kryoSerializer.newInstance();
        final ClassTag<T> tag = ClassTag$.MODULE$.apply(inputClazz);
        return sparkSerializer.deserialize(sparkSerializer.serialize(input, tag), tag);
    }
}
