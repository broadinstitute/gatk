package org.broadinstitute.hellbender.testutils;


import org.apache.spark.SparkConf;
import org.apache.spark.serializer.KryoSerializer;
import org.apache.spark.serializer.SerializerInstance;
import org.broadinstitute.hellbender.utils.Utils;
import scala.reflect.ClassTag;
import scala.reflect.ClassTag$;

import java.io.*;

public final class SparkTestUtils {
    private SparkTestUtils() {}

    /**
     * Takes an input object and returns the value of the object after it has been serialized and then deserialized in Kryo.
     * Requires the class of the input object as a parameter because it's not generally possible to get the class of a
     * generified method parameter with reflection.
     *
     * @param input instance of inputClazz.  Never {@code null}
     * @param inputClazz class to cast input
     * @param conf Spark configuration to test
     * @param <T> class to attempt.  Same or subclass of inputClazz
     * @return serialized and deserialized instance of input.  Throws exception if serialization round trip fails.
     */
    public static <T> T roundTripInKryo(final T input, final Class<?> inputClazz, final SparkConf conf) {
        Utils.nonNull(input);
        final KryoSerializer kryoSerializer = new KryoSerializer(conf);
        final SerializerInstance sparkSerializer = kryoSerializer.newInstance();
        final ClassTag<T> tag = ClassTag$.MODULE$.apply(inputClazz);
        return sparkSerializer.deserialize(sparkSerializer.serialize(input, tag), tag);
    }

    /**
     * Takes an input object and returns the value of the object after it has been serialized and then deserialized
     * using Java's built in serialization.
     *
     * @param input an object to be serialized.  Never {@code null}
     * @return serialized and deserialized instance of input, may throw if serialization fails
     */
    @SuppressWarnings("unchecked")
    public static <T> T roundTripThroughJavaSerialization(T input) throws IOException, ClassNotFoundException {
        Utils.nonNull(input);
        final byte[] serializedBytes;
        try(final ByteArrayOutputStream baos = new ByteArrayOutputStream();
            final ObjectOutputStream out = new ObjectOutputStream(baos)){
            out.writeObject(input);
            serializedBytes = baos.toByteArray();
        }

        try (final ByteArrayInputStream bais = new ByteArrayInputStream(serializedBytes);
             final ObjectInputStream in = new ObjectInputStream(bais)){

            return (T) in.readObject();
        }
    }
}
