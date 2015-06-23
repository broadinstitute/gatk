package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.api.client.util.Sets;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.values.*;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.util.*;


/**
 * RemoveDuplicateReadVariantPairs removes duplicates <GATKRead,Variant> pairs, then groups by GATKRead to produce
 * <GATKRead,Iterable<Variant>>. Conceptually, this is just
 * input.apply(RemoveDuplicates.<KV<GATKRead, Variant>>create()).apply(GroupByKey.<GATKRead, Variant>create());
 *
 * Unfortunately, GATKRead cannot be deterministically encoded, so we use UUID as surrogates for GATKRead and Variant
 * to remove duplicate KV pair values. Then, we have the laborious process of adding the "real" data back in.
 *
 * This currently requires 9 transforms.
 * (1) pTransform to get PCollection<KV<UUID, UUID>>
 * (2) remove dupes of (1)
 * (3) GroupByKey of (2) (produces PCollection<UUID,Iterable<UUID>>)
 * (4) create PCollection<UUID, KV<UUID, Variant>>
 * (5) create KeyedPCollectionTuple of (3) and (4)
 * (6) use (5) to create PCollection<UUID, Iterable<Variant>>
 * (7) use  pKvAB to create PCollection<KV<UUID, GATKRead>>
 * (8) create KeyedPCollectionTuple of (6) and (7)
 * (9) join Alpha back in (using (8)) to get PCollection<GATKRead, Iterable<Variant>>
 *
 */
public class RemoveDuplicateReadVariantPairs extends PTransform<PCollection<KV<GATKRead,Variant>>, PCollection<KV<GATKRead, Iterable<Variant>>>> {
    private static final long serialVersionUID = 1L;

    @Override
    public PCollection<KV<GATKRead, Iterable<Variant>>> apply(PCollection<KV<GATKRead, Variant>> input) {
        PCollection<KV<UUID, UUID>> readVariantUuids = input.apply(ParDo.of(new StripToUUIDs())).setName("RemoveDuplicateReadVariantPairs_stripToUUIDs");

        // RemoveDuplicates removes duplicate KV pairs (not just the keys).
        PCollection<KV<UUID, Iterable<UUID>>> readsIterableVariants = readVariantUuids.apply(RemoveDuplicates.<KV<UUID, UUID>>create()).apply(GroupByKey.<UUID, UUID>create());

        // Now join the stuff we care about back in.
        // Start by making KV<ReadUUID, KV<VariantUUID, Variant>> and joining that to KV<ReadUUID, Iterable<VariantUUID>>
        PCollection<KV<UUID, Iterable<Variant>>> matchedVariants =
                RemoveReadVariantDupesUtility.addBackVariants(input, readsIterableVariants);

        // And now, we do the same song and dance to get the Reads back in.
        return RemoveReadVariantDupesUtility.addBackReads(input, matchedVariants);
    }

    static class StripToUUIDs extends DoFn<KV<GATKRead, Variant>, KV<UUID, UUID>> {
        private static final long serialVersionUID = 1L;
        @Override
        public void processElement(ProcessContext c) throws Exception {
            c.output(KV.of(c.element().getKey().getUUID(), c.element().getValue().getUUID()));
        }
    }

    static class RemoveReadVariantDupesUtility {
        public static PCollection<KV<UUID, Iterable<Variant>>> addBackVariants(PCollection<KV<GATKRead, Variant>> readVariants, PCollection<KV<UUID, Iterable<UUID>>> readsUUIDsIterableVariantUUIDs) {
            // Now join the stuff we care about back in.
            // Start by making KV<ReadUUID, KV<VariantUUID, Variant>> and joining that to KV<ReadUUID, Iterable<VariantUUID>>
            PCollection<KV<UUID, KV<UUID, Variant>>> variantToBeJoined =
                    readVariants.apply(ParDo.of(new DoFn<KV<GATKRead, Variant>, KV<UUID, KV<UUID, Variant>>>() {
                        private static final long serialVersionUID = 1L;
                        @Override
                        public void processElement(ProcessContext c) throws Exception {
                            c.output(KV.of(c.element().getKey().getUUID(), KV.of(c.element().getValue().getUUID(), c.element().getValue())));
                        }
                    })).setName("RemoveReadVariantDupesUtility_VariantToBeJoined");

            // Another coGroupBy...
            final TupleTag<Iterable<UUID>> iterableUuidTag = new TupleTag<>(); // Iterable<VariantUUID>
            final TupleTag<KV<UUID, Variant>> variantUuidTag = new TupleTag<>(); // KV<VariantUUID, Variant>

            PCollection<KV<UUID, CoGbkResult>> coGbkAgain = KeyedPCollectionTuple
                    .of(iterableUuidTag, readsUUIDsIterableVariantUUIDs)
                    .and(variantUuidTag, variantToBeJoined).apply(CoGroupByKey.<UUID>create());

            // For each ReadUUID, get all of the variants that have a VariantUUID
            return coGbkAgain.apply(ParDo.of(
                    new DoFn<KV<UUID, CoGbkResult>, KV<UUID, Iterable<Variant>>>() {
                        private static final long serialVersionUID = 1L;

                        @Override
                        public void processElement(ProcessContext c) throws Exception {
                            // kVariants UUID is the UUID for the Variant
                            Iterable<KV<UUID, Variant>> kVariants = c.element().getValue().getAll(variantUuidTag);
                            // kUUIDss is two deep iterable of VariantUUIDs
                            Iterable<Iterable<UUID>> kUUIDss = c.element().getValue().getAll(iterableUuidTag);
                            // For every UUID that's left, keep the variant.
                            Set<UUID> uuidHashSet = Sets.newHashSet();
                            for (Iterable<UUID> uuids : kUUIDss) {
                                for (UUID uuid : uuids) {
                                    uuidHashSet.add(uuid);
                                }
                            }

                            // Now find every variant for each UUID.
                            Map<UUID, Variant> variantMap = Maps.newHashMap();
                            for (KV<UUID,Variant> uVariants : kVariants) {
                                if (uuidHashSet.contains(uVariants.getKey()) &&
                                        !variantMap.containsKey(uVariants.getKey())) {
                                    variantMap.put(uVariants.getKey(), uVariants.getValue());
                                }
                            }
                            List<Variant> iVariants = Lists.newArrayList();
                            Iterator<Map.Entry<UUID, Variant>> iterator = variantMap.entrySet().iterator();
                            while (iterator.hasNext()) {
                                iVariants.add(iterator.next().getValue());
                                iterator.remove();
                            }

                            c.output(KV.of(c.element().getKey(), iVariants));
                        }
                    })).setName("RemoveReadVariantDupesUtility_CoGroupBy");

        }

        public static PCollection<KV<GATKRead, Iterable<Variant>>> addBackReads(PCollection<KV<GATKRead, Variant>> readVariants, PCollection<KV<UUID, Iterable<Variant>>> matchedVariants) {
            // And now, we do the same song and dance to get the Reads back in.
            final TupleTag<GATKRead> justReadTag = new TupleTag<>();
            final TupleTag<Iterable<Variant>> iterableVariant = new TupleTag<>();

            PCollection<KV<UUID, GATKRead>> kReads = readVariants.apply(Keys.<GATKRead>create()).apply(new KeyReadsByUUID());
            PCollection<KV<UUID, CoGbkResult>> coGbkLast = KeyedPCollectionTuple
                    .of(justReadTag, kReads)
                    .and(iterableVariant, matchedVariants).apply(CoGroupByKey.<UUID>create());

            return coGbkLast.apply(ParDo.of(new DoFn<KV<UUID, CoGbkResult>, KV<GATKRead, Iterable<Variant>>>() {
                private static final long serialVersionUID = 1L;
                @Override
                public void processElement(ProcessContext c) throws Exception {
                    Iterable<GATKRead> iReads = c.element().getValue().getAll(justReadTag);
                    // We only care about the first read (the rest are the same.
                    Iterable<Iterable<Variant>> variants = c.element().getValue().getAll(iterableVariant);
                    List<GATKRead> reads = Lists.newArrayList();
                    for (GATKRead r : iReads) {
                        reads.add(r);
                    }
                    if (reads.isEmpty()) {
                        throw new GATKException("no reads found");
                    }

                    for (Iterable<Variant> v : variants) {
                        c.output(KV.of(reads.get(0), v));
                    }
                }
            })).setName("RemoveDuplicatePairedReadVariants_addBackReads");
        }
    }
}
