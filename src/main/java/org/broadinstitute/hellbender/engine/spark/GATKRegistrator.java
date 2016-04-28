package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.google.api.services.genomics.model.Read;
import com.esotericsoftware.kryo.serializers.FieldSerializer;
import de.javakaffee.kryoserializers.UnmodifiableCollectionsSerializer;
import htsjdk.samtools.SAMRecord;
import org.apache.spark.serializer.KryoRegistrator;
import org.bdgenomics.adam.serialization.ADAMKryoRegistrator;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.PairedEnds;

import java.util.Collections;

/**
 * GATKRegistrator registers Serializers for our project. We need a JsonSerializer for the Google Genomics classes
 * and UnmodifiableCollectionsSerializer from a bug in the version of Kryo we're on.
 */
public class GATKRegistrator implements KryoRegistrator {
    private static KryoRegistrator wrappedRegistrator = new DefaultRegistrator();

    public static KryoRegistrator getCurrentRegistrator() { return wrappedRegistrator; }

    public static KryoRegistrator setCurrentRegistrator( KryoRegistrator newRegistrator ) {
        final KryoRegistrator oldRegistrator = wrappedRegistrator;
        wrappedRegistrator = newRegistrator != null ? newRegistrator : new DefaultRegistrator();
        return oldRegistrator;
    }

    @Override
    public void registerClasses( Kryo kryo ) { wrappedRegistrator.registerClasses(kryo); }

    public static final class DefaultRegistrator implements KryoRegistrator {
        private static final ADAMKryoRegistrator adamRegistrator = new ADAMKryoRegistrator();

        @Override
        public void registerClasses(Kryo kryo) {

            // JsonSerializer is needed for the Google Genomics classes like Read and Reference.
            kryo.register(Read.class, new JsonSerializer<Read>());
            // htsjdk.variant.variantcontext.CommonInfo has a Map<String, Object> that defaults to
            // a Collections.unmodifiableMap. This can't be handled by the version of kryo used in Spark, it's fixed
            // in newer versions (3.0.x), but we can't use those because of incompatibility with Spark. We just include the
            // fix here.
            // We are tracking this issue with (#874)
            kryo.register(Collections.unmodifiableMap(Collections.emptyMap()).getClass(), new UnmodifiableCollectionsSerializer());

            kryo.register(Collections.unmodifiableList(Collections.emptyList()).getClass(), new UnmodifiableCollectionsSerializer());

            kryo.register(SAMRecordToGATKReadAdapter.class, new SAMRecordToGATKReadAdapterSerializer());

            kryo.register(SAMRecord.class, new SAMRecordSerializer());

            //register to avoid writing the full name of this class over and over
            kryo.register(PairedEnds.class, new FieldSerializer<>(kryo, PairedEnds.class));

            // register the ADAM data types using Avro serialization, including:
            //     AlignmentRecord
            //     Genotype
            //     Variant
            //     DatabaseVariantAnnotation
            //     NucleotideContigFragment
            //     Contig
            //     StructuralVariant
            //     VariantCallingAnnotations
            //     VariantEffect
            //     DatabaseVariantAnnotation
            //     Dbxref
            //     Feature
            //     ReferencePosition
            //     ReferencePositionPair
            //     SingleReadBucket
            //     IndelRealignmentTarget
            //     TargetSet
            //     ZippedTargetSet
            adamRegistrator.registerClasses(kryo);
        }
    }
}
