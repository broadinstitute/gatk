package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.serializer.SerializerInstance;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

public class SparkTestUtils {

    /**
     * Takes an input object and returns the value of the object after it has been serialized and then deserialized in Kryo.
     * Requires the class of the input object as a parameter because it's not gererally possible to get the class of a
     * generified method parameter with reflection.
     * @param input The object to be serialized and deserialized.
     * @param inputClazz The class of the input object.
     * @param conf SparkConf with any necessary Kryo registrators defined
     * @return
     */
    public static <T> T roundTripInKryo(T input, Class<?> inputClazz, final SparkConf conf) {
        KryoSerializer kryoSerializer = new KryoSerializer(conf);
        SerializerInstance sparkSerializer = kryoSerializer.newInstance();
        final ClassTag<T> tag = ClassTag$.MODULE$.apply(inputClazz);
        return sparkSerializer.deserialize(sparkSerializer.serialize(input, tag), tag);
    }
}
