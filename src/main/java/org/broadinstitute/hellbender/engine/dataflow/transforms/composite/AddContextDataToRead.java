package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.*;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.dev.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPIMetadata;
import org.broadinstitute.hellbender.engine.dataflow.datasources.RefAPISource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantsDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.transforms.KeyReadsByUUID;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.List;
import java.util.UUID;


/**
 * AddContextDataToRead pairs reference bases and overlapping variants with each GATKRead in the PCollection input.
 * The variants are obtained from a local file (later a GCS Bucket). The reference bases come from the Google Genomics API.
 *
 * This transform is intended for direct use in pipelines.
 */
public class AddContextDataToRead {
    public static PCollection<KV<GATKRead, ReadContextData>> add(PCollection<GATKRead> pReads, RefAPIMetadata refAPIMetadata, VariantsDataflowSource variantsDataflowSource) {
        PCollection<Variant> pVariants = variantsDataflowSource.getAllVariants();
        PCollection<KV<GATKRead, Iterable<Variant>>> kvReadVariants = KeyVariantsByRead.key(pVariants, pReads);
        PCollection<KV<GATKRead, ReferenceBases>> kvReadRefBases =
                RefBasesForReads.addBases(pReads, refAPIMetadata);
        return join(pReads, kvReadRefBases, kvReadVariants);

    }

    /**
     * Preconditions per read: Exactly one ReferenceBases object and zero or one iterable of Variants.
     */
    protected static PCollection<KV<GATKRead, ReadContextData>> join(PCollection<GATKRead> pReads, PCollection<KV<GATKRead, ReferenceBases>> kvReadRefBases, PCollection<KV<GATKRead, Iterable<Variant>>> kvReadVariants) {
        boolean assertsEnabled = false;
        assert assertsEnabled = true; // Intentional side-effect!!!
        // Now assertsEnabled is set to the correct value
        if (assertsEnabled) {
            assertSameReads(pReads, kvReadRefBases, kvReadVariants);
        }
        // We could add a check that all of the reads in kvReadRefBases, pVariants, and pReads are the same.
        PCollection<KV<UUID, Iterable<Variant>>> UUIDVariants = kvReadVariants.apply(ParDo.named("KvUUIDVariants").of(new DoFnWLog<KV<GATKRead, Iterable<Variant>>, KV<UUID, Iterable<Variant>>>("KvUUIDVariants") {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                GATKRead r = c.element().getKey();
                Iterable<Variant> variants = c.element().getValue();
                c.output(KV.of(r.getUUID(), variants));
            }
        })).setName("KvUUIDiVariants");

        PCollection<KV<UUID, ReferenceBases>> UUIDRefBases = kvReadRefBases.apply(ParDo.named("KvUUIDRefBases").of(new DoFnWLog<KV<GATKRead, ReferenceBases>, KV<UUID, ReferenceBases>>("KvUUIDRefBases") {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                GATKRead r = c.element().getKey();
                c.output(KV.of(r.getUUID(), c.element().getValue()));
            }
        })).setName("KvUUIDRefBases");

        final TupleTag<Iterable<Variant>> variantTag = new TupleTag<>();
        final TupleTag<ReferenceBases> referenceTag = new TupleTag<>();

        PCollection<KV<UUID, CoGbkResult>> coGbkInput = KeyedPCollectionTuple
                .of(variantTag, UUIDVariants)
                .and(referenceTag, UUIDRefBases).apply(CoGroupByKey.<UUID>create());

        PCollection<KV<UUID, ReadContextData>> UUIDcontext = coGbkInput.apply(ParDo.named("kVUUIDReadContextData").of(new DoFnWLog<KV<UUID, CoGbkResult>, KV<UUID, ReadContextData>>("kVUUIDReadContextData") {
            private static final long serialVersionUID = 1L;

            @Override
            public void processElement(ProcessContext c) throws Exception {
                Iterable<Iterable<Variant>> variants = c.element().getValue().getAll(variantTag);
                Iterable<ReferenceBases> referenceBases = c.element().getValue().getAll(referenceTag);

                List<Iterable<Variant>> vList = Lists.newArrayList(variants);
                if (vList.isEmpty()) {
                    vList.add(Lists.newArrayList());
                }
                if (vList.size() > 1) {
                    throw new GATKException("expected one list of variants, got " + vList.size());
                }
                List<ReferenceBases> bList = Lists.newArrayList(referenceBases);
                if (bList.size() != 1) {
                    throw new GATKException("expected one ReferenceBases, got " + bList.size());
                }
                c.output(KV.of(c.element().getKey(), new ReadContextData(bList.get(0), vList.get(0))));
            }
        })).setName("kVUUIDReadContextData");

        // Now add the reads back in: key by UUID, then join back in.
        PCollection<KV<UUID, GATKRead>> UUIDRead = pReads.apply(new KeyReadsByUUID());
        final TupleTag<GATKRead> readTag = new TupleTag<>();
        final TupleTag<ReadContextData> contextDataTag = new TupleTag<>();

        PCollection<KV<UUID, CoGbkResult>> coGbkfull = KeyedPCollectionTuple
                .of(readTag, UUIDRead)
                .and(contextDataTag, UUIDcontext).apply(CoGroupByKey.<UUID>create());

        return coGbkfull.apply(ParDo.named("joinReadwithContextData").of(new DoFnWLog<KV<UUID, CoGbkResult>, KV<GATKRead, ReadContextData>>("joinReadwithContextData") {
            private static final long serialVersionUID = 1L;

            @Override
            public void processElement(ProcessContext c) throws Exception {
                Iterable<GATKRead> reads = c.element().getValue().getAll(readTag);
                Iterable<ReadContextData> contextDatas = c.element().getValue().getAll(contextDataTag);

                List<GATKRead> rList = Lists.newArrayList();
                for (GATKRead r : reads) {
                    rList.add(r);
                }
                if (rList.size() != 1) {
                    throw new GATKException("expected one Read, got " + rList.size());
                }
                List<ReadContextData> cList = Lists.newArrayList();
                for (ReadContextData cd : contextDatas) {
                    cList.add(cd);
                }
                if (cList.size() != 1) {
                    throw new GATKException("expected one ReadContextData, got " + cList.size());
                }

                c.output(KV.of(rList.get(0), cList.get(0)));
            }
        })).setName("joinReadwithContextData");
    }

    static void assertSameReads(PCollection<GATKRead> pReads, PCollection<KV<GATKRead, ReferenceBases>> kvReadRefBases, PCollection<KV<GATKRead, Iterable<Variant>>> kvReadVariants) {
        // Because dataflow is inherently parallel, determining overlap is a little hard. Here, we
        // (1) Get the UUIDs for each pCollection and count the UUIDs
        // (2) Put the three collections into a single collection, find, and count the unique UUIDs
        // (3) Verify that all of the counts match.
        PCollection<UUID> readUUIDs = pReads.apply(new StripToUUID());
        PCollection<UUID> readRefBasesUUIDs = kvReadRefBases.apply(Keys.<GATKRead>create()).apply(new StripToUUID());
        PCollection<UUID> readVariantsUUIDs = kvReadVariants.apply(Keys.<GATKRead>create()).apply(new StripToUUID());
        PCollectionView<Long> readCount = readUUIDs.apply(Count.<UUID>globally()).apply(View.asSingleton());
        PCollectionView<Long> readRefBasesCount = readUUIDs.apply(Count.<UUID>globally()).apply(View.asSingleton());
        PCollectionView<Long> readVariantsCount = readUUIDs.apply(Count.<UUID>globally()).apply(View.asSingleton());

        PCollectionList<UUID> allUUIDs = PCollectionList.of(readUUIDs).and(readRefBasesUUIDs).and(readVariantsUUIDs);
        PCollection<UUID> flatUUIDs = allUUIDs.apply(Flatten.<UUID>pCollections());
        PCollection<UUID> uniqueUUIDs = flatUUIDs.apply(RemoveDuplicates.<UUID>create());
        PCollection<Long> uniqueUUIDCount = uniqueUUIDs.apply(Count.<UUID>globally());

        uniqueUUIDCount.apply(ParDo.named("assertSameReads")
            .withSideInputs(readCount, readRefBasesCount, readVariantsCount)
            .of(new DoFnWLog<Long, Long>("assertSameReads") {
            private static final long serialVersionUID = 1L;

            @Override
            public void processElement(ProcessContext c) throws Exception {
                Long reads = c.sideInput(readCount);
                Long readRefBases = c.sideInput(readRefBasesCount);
                Long readVariants = c.sideInput(readVariantsCount);
                assert c.element().equals(reads);
                assert c.element().equals(readRefBases);
                assert c.element().equals(readVariants);
            }
        }));
    }

    static class StripToUUID extends PTransform<PCollection<GATKRead>, PCollection<UUID>> {
        private static final long serialVersionUID = 1L;

        @Override
        public PCollection<UUID> apply(PCollection<GATKRead> input) {
            return input.apply(ParDo.named("StripToUUID").of(new DoFnWLog<GATKRead, UUID>("StripToUUID") {
                private static final long serialVersionUID = 1L;

                @Override
                public void processElement(ProcessContext c) throws Exception {
                    c.output(c.element().getUUID());
                }
            }));
        }
    }
}