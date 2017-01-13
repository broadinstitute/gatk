package org.broadinstitute.hellbender.utils.spark;

import com.esotericsoftware.kryo.Kryo;
import de.javakaffee.kryoserializers.UnmodifiableCollectionsSerializer;
import org.apache.spark.serializer.KryoRegistrator;

/**
 * A Kryo registrator for unmodifiable collections serializers
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class UnmodifiableCollectionsRegistrator implements KryoRegistrator {
    @Override
    public void registerClasses(Kryo kryo) {
        UnmodifiableCollectionsSerializer.registerSerializers(kryo);
    }
}
