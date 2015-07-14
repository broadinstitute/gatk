package org.broadinstitute.hellbender.engine.dataflow.transforms.composite;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.TupleTag;
import com.google.common.collect.Lists;
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
    public static PCollection<KV<GATKRead, ReadContextData>> add(PCollection<GATKRead> pReads, RefAPISource refAPISource, RefAPIMetadata refAPIMetadata, VariantsDataflowSource variantsDataflowSource) {
        PCollection<Variant> pVariants = variantsDataflowSource.getAllVariants();
        PCollection<KV<GATKRead, Iterable<Variant>>> kvReadVariants = KeyVariantsByRead.key(pVariants, pReads);
        PCollection<KV<GATKRead, ReferenceBases>> kvReadRefBases =
                RefBasesForReads.addBases(pReads, refAPISource, refAPIMetadata);
        return join(pReads, kvReadRefBases, kvReadVariants);

    }

    /**
     * Preconditions per read: Exactly one ReferenceBases object and zero or one iterable of Variants.
     */
    protected static PCollection<KV<GATKRead, ReadContextData>> join(PCollection<GATKRead> pReads, PCollection<KV<GATKRead, ReferenceBases>> kvReadRefBases, PCollection<KV<GATKRead, Iterable<Variant>>> kvReadVariants) {
        // We could add a check that all of the reads in kvReadRefBases, pVariants, and pReads are the same.
        PCollection<KV<UUID, Iterable<Variant>>> UUIDVariants = kvReadVariants.apply(ParDo.of(new DoFn<KV<GATKRead, Iterable<Variant>>, KV<UUID, Iterable<Variant>>>() {
            private static final long serialVersionUID = 1L;
            @Override
            public void processElement(ProcessContext c) throws Exception {
                GATKRead r = c.element().getKey();
                Iterable<Variant> variants = c.element().getValue();
                c.output(KV.of(r.getUUID(), variants));
            }
        })).setName("KvUUIDiVariants");

        PCollection<KV<UUID, ReferenceBases>> UUIDRefBases = kvReadRefBases.apply(ParDo.of(new DoFn<KV<GATKRead, ReferenceBases>, KV<UUID, ReferenceBases>>() {
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

        PCollection<KV<UUID, ReadContextData>> UUIDcontext = coGbkInput.apply(ParDo.of(new DoFn<KV<UUID, CoGbkResult>, KV<UUID, ReadContextData>>() {
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

        return coGbkfull.apply(ParDo.of(new DoFn<KV<UUID, CoGbkResult>, KV<GATKRead, ReadContextData>>() {
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

}