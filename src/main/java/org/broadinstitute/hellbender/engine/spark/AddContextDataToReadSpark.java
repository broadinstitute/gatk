package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.util.List;
import java.util.Map;

/**
 * AddContextDataToRead pairs reference bases and overlapping variants with each GATKRead in the PCollection input.
 * The variants are obtained from a local file (later a GCS Bucket). The reference bases come from the Google Genomics API.
 *
 * This transform is intended for direct use in pipelines.
 *
 * The reference bases paired with each read can be customized by passing in a reference window function
 * inside the {@link ReferenceDataflowSource} argument to {@link #add}. See
 * {@link org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions} for examples.
 */
public class AddContextDataToReadSpark {
    public static JavaPairRDD<GATKRead, ReadContextData> add(
            final JavaRDD<GATKRead> reads, final ReferenceDataflowSource referenceDataflowSource,
            final JavaRDD<Variant> variants) {
        // Join Reads and Variants, Reads and ReferenceBases
        JavaPairRDD<GATKRead, Iterable<Variant>> readiVariants = JoinReadsWithVariants.join(reads, variants);
        JavaPairRDD<GATKRead, ReferenceBases> readRefBases = JoinReadsWithRefBases.addBases(referenceDataflowSource, reads);

        // For testing we want to know that the reads from the KVs coming back from JoinReadsWithVariants.Join
        // and JoinReadsWithRefBases.Pair are the same reads from "reads".
        boolean assertsEnabled = false;
        assert assertsEnabled = true; // Intentional side-effect!!!
        // Now assertsEnabled is set to the correct value
        if (assertsEnabled) {
            assertSameReads(reads, readRefBases, readiVariants);
        }

        JavaPairRDD<GATKRead, Tuple2<Iterable<Iterable<Variant>>, Iterable<ReferenceBases>>> cogroup = readiVariants.cogroup(readRefBases);
        return cogroup.mapToPair(in -> {
            List<Iterable<Variant>> liVariants = Lists.newArrayList(in._2()._1());
            List<Variant> lVariants = Lists.newArrayList();
            if (!liVariants.isEmpty()) {
                final Iterable<Variant> iVariant = Iterables.getOnlyElement(in._2()._1());
                // It's possible for the iVariant to contain only a null variant, we don't
                // want to add that to the ReadContextData.
                final Variant next = iVariant.iterator().next();
                if (next != null) {
                    lVariants = Lists.newArrayList(iVariant);
                }
            }

            ReferenceBases refBases = Iterables.getOnlyElement(in._2()._2());
            ReadContextData readContextData = new ReadContextData(refBases, lVariants);
            return new Tuple2<>(in._1(), readContextData);
        });
    }

    private static void assertSameReads(final JavaRDD<GATKRead> reads,
                                        final JavaPairRDD<GATKRead, ReferenceBases> readRefBases,
                                        final JavaPairRDD<GATKRead, Iterable<Variant>> readiVariants) {
        List<GATKRead> vReads = reads.collect();
        Map<GATKRead, ReferenceBases> vReadRef = readRefBases.collectAsMap();
        Map<GATKRead, Iterable<Variant>> vReadVariant = readiVariants.collectAsMap();

        // This assumes all rdds doesn't have any duplicates.
        JavaRDD<GATKRead> refBasesReads = readRefBases.keys();
        JavaRDD<GATKRead> variantsReads = readiVariants.keys();
        JavaRDD<GATKRead> distinctReads = reads.intersection(refBasesReads).intersection(variantsReads);

        long counts = reads.count();
        assert counts == distinctReads.count();
        assert counts == refBasesReads.count();
        assert counts == variantsReads.count();

        assert vReadRef.size() == vReads.size();
        assert vReadVariant.size() == vReads.size();
    }

}

