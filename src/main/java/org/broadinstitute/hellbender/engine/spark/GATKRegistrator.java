package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.serializers.FieldSerializer;
import com.esotericsoftware.kryo.serializers.JavaSerializer;
import com.google.api.services.genomics.model.Read;
import de.javakaffee.kryoserializers.UnmodifiableCollectionsSerializer;
import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.FastGenotype;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
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

    private ADAMKryoRegistrator ADAMregistrator;

    public GATKRegistrator() {
        this.ADAMregistrator = new ADAMKryoRegistrator();
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    @Override
    public void registerClasses(Kryo kryo) {

        // JsonSerializer is needed for the Google Genomics classes like Read and Reference.
        kryo.register(Read.class, new JsonSerializer<Read>());

        //relatively inefficient serialization of Collections created with Collections.nCopies(), without this
        //any Collection created with Collections.nCopies fails to serialize at run time
        kryo.register(Collections.nCopies(2, "").getClass(), new FieldSerializer<>(kryo, Collections.nCopies(2, "").getClass()));

        // htsjdk.variant.variantcontext.CommonInfo has a Map<String, Object> that defaults to
        // a Collections.unmodifiableMap. This can't be handled by the version of kryo used in Spark, it's fixed
        // in newer versions (3.0.x), but we can't use those because of incompatibility with Spark. We just include the
        // fix here.
        // We are tracking this issue with (#874)
        kryo.register(Collections.unmodifiableMap(Collections.EMPTY_MAP).getClass(), new UnmodifiableCollectionsSerializer());

        kryo.register(Collections.unmodifiableList(Collections.EMPTY_LIST).getClass(), new UnmodifiableCollectionsSerializer());

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
        ADAMregistrator.registerClasses(kryo);
    }
}
