package org.broadinstitute.hellbender.engine.spark;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.serializers.FieldSerializer;
import com.google.common.collect.ImmutableMap;
import de.javakaffee.kryoserializers.UnmodifiableCollectionsSerializer;
import de.javakaffee.kryoserializers.guava.ImmutableMapSerializer;
import htsjdk.samtools.*;
import org.apache.spark.serializer.KryoRegistrator;
import org.bdgenomics.adam.serialization.ADAMKryoRegistrator;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords.*;

import java.util.Collections;

/**
 * GATKRegistrator registers Serializers for our project. We need a JsonSerializer for the Google Genomics classes
 * and UnmodifiableCollectionsSerializer from a bug in the version of Kryo we're on.
 */
public class GATKRegistrator implements KryoRegistrator {

    private final ADAMKryoRegistrator ADAMregistrator = new ADAMKryoRegistrator();

    public GATKRegistrator() {}

    @Override
    public void registerClasses(Kryo kryo) {

        registerGATKClasses(kryo);


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

        //do this before and after ADAM to try and force our registrations to win out
        registerGATKClasses(kryo);
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private void registerGATKClasses(Kryo kryo) {
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

        kryo.register(ImmutableMap.of().getClass(), new ImmutableMapSerializer());
        kryo.register(ImmutableMap.of("one","element").getClass(), new ImmutableMapSerializer());
        kryo.register(ImmutableMap.of("map","with","multiple","elements").getClass(), new ImmutableMapSerializer());

        kryo.register(SAMRecordToGATKReadAdapter.class, new SAMRecordToGATKReadAdapterSerializer());

        kryo.register(SAMRecord.class, new SAMRecordSerializer());
        kryo.register(BAMRecord.class, new SAMRecordSerializer());

        kryo.register(SAMFileHeader.class);
        kryo.register(SAMFileHeader.GroupOrder.class);
        kryo.register(SAMFileHeader.SortOrder.class);
        kryo.register(SAMProgramRecord.class);
        kryo.register(SAMReadGroupRecord.class);
        kryo.register(EmptyFragment.class, new FieldSerializer(kryo, EmptyFragment.class));
        kryo.register(Fragment.class, new FieldSerializer(kryo, Fragment.class));
        kryo.register(Pair.class, new Pair.Serializer());
        kryo.register(Passthrough.class, new FieldSerializer(kryo, Passthrough.class));
        kryo.register(MarkDuplicatesSparkUtils.IndexPair.class, new FieldSerializer(kryo, MarkDuplicatesSparkUtils.IndexPair.class));
        kryo.register(ReadsKey.class, new FieldSerializer(kryo, ReadsKey.class));
        kryo.register(ReadsKey.KeyForFragment.class, new FieldSerializer(kryo, ReadsKey.KeyForFragment.class));
        kryo.register(ReadsKey.KeyForPair.class, new FieldSerializer(kryo, ReadsKey.KeyForPair.class));
    }
}
