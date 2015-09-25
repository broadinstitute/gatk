package org.broadinstitute.hellbender.engine.spark;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.util.List;
import java.util.NoSuchElementException;

/**
 * AddContextDataToRead pairs reference bases and overlapping variants with each GATKRead in the RDD input.
 * The variants are obtained from a local file (later a GCS Bucket). The reference bases come from the Google Genomics API.
 *
 * This transform is intended for direct use in pipelines.
 *
 * This transform will filter out any unmapped reads.
 *
 * The reference bases paired with each read can be customized by passing in a reference window function
 * inside the {@link ReferenceDataflowSource} argument to {@link #add}. See
 * {@link org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions} for examples.
 */
public class AddContextDataToReadSpark {
    public static JavaPairRDD<GATKRead, ReadContextData> add(
            final JavaRDD<GATKRead> reads, final ReferenceDataflowSource referenceDataflowSource,
            final JavaRDD<Variant> variants) {

        // This transform can currently only handle mapped reads
        JavaRDD<GATKRead> mappedReads = reads.filter(read -> ReadFilterLibrary.MAPPED.test(read));

        // Join Reads and Variants, Reads and ReferenceBases
        JavaPairRDD<GATKRead, Iterable<Variant>> readiVariants = JoinReadsWithVariants.join(mappedReads, variants);
        JavaPairRDD<GATKRead, ReferenceBases> readRefBases = JoinReadsWithRefBases.addBases(referenceDataflowSource, mappedReads);

        // For testing we want to know that the reads from the KVs coming back from JoinReadsWithVariants.Join
        // and JoinReadsWithRefBases.Pair are the same reads from "reads".
        boolean assertsEnabled = false;
        assert assertsEnabled = true; // Intentional side-effect!!!
        // Now assertsEnabled is set to the correct value
        if (assertsEnabled) {
            assertSameReads(mappedReads, readRefBases, readiVariants);
        }

        JavaPairRDD<GATKRead, Tuple2<Iterable<Iterable<Variant>>, Iterable<ReferenceBases>>> cogroup = readiVariants.cogroup(readRefBases);
        return cogroup.mapToPair(in -> {
            ReadContextData readContextData = null;
            try {
                List<Variant> lVariants = makeListFromIterableIterable(in._2()._1());

                ReferenceBases refBases = Iterables.getOnlyElement(in._2()._2());
                readContextData = new ReadContextData(refBases, lVariants);
            } catch(NoSuchElementException e) {
                throw new GATKException.ShouldNeverReachHereException(e);
            }
            return new Tuple2<>(in._1(), readContextData);
        });
    }

    private static <T> List<T> makeListFromIterableIterable(Iterable<Iterable<T>> iterables) {
        List<Iterable<T>> listIterableT = Lists.newArrayList(iterables);
        List<T> listT = Lists.newArrayList();
        if (!listIterableT.isEmpty()) {
            final Iterable<T> iterableT = Iterables.getOnlyElement(iterables);
            // It's possible for the iterableT to contain only a null T, we don't
            // want to include that.
            final T next = iterableT.iterator().next();
            if (next != null) {
                listT = Lists.newArrayList(iterableT);
            }
        }
        return listT;

    }

    private static void assertSameReads(final JavaRDD<GATKRead> reads,
                                        final JavaPairRDD<GATKRead, ReferenceBases> readRefBases,
                                        final JavaPairRDD<GATKRead, Iterable<Variant>> readiVariants) {

        // We want to verify that the reads are the same for each collection and that there are no duplicates
        // in any collection.

        // Collect all reads (with potential duplicates) in allReads. We expect there to be 3x the unique reads.
        // Verify that there are 3x the distinct reads and the reads count for each collection match
        // the distinct reads count.
        // We should also check that the reference bases and variants are correctly paired with the reads. See
        // issue (#873).
        JavaRDD<GATKRead> refBasesReads = readRefBases.keys();
        JavaRDD<GATKRead> variantsReads = readiVariants.keys();
        JavaRDD<GATKRead> allReads = reads.union(refBasesReads).union(variantsReads);
        long allReadsCount = allReads.count();
        long distinctReads = allReads.distinct().count();

        assert 3*distinctReads == allReadsCount;
        assert distinctReads == reads.count();
        assert distinctReads == refBasesReads.count();
        assert distinctReads == variantsReads.count();
    }

}

