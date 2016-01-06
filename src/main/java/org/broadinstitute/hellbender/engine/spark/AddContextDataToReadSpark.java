package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

/**
 * AddContextDataToRead pairs reference bases and overlapping variants with each GATKRead in the RDD input.
 * The variants are obtained from a local file (later a GCS Bucket). The reference bases come from the Google Genomics API.
 *
 * This transform is intended for direct use in pipelines.
 *
 * This transform will filter out any unmapped reads.
 *
 * The reference bases paired with each read can be customized by passing in a reference window function
 * inside the {@link org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource} argument to {@link #add}. See
 * {@link org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions} for examples.
 */
public class AddContextDataToReadSpark {
    public static JavaPairRDD<GATKRead, ReadContextData> add(
            final JavaRDD<GATKRead> reads, final ReferenceMultiSource referenceDataflowSource,
            final JavaRDD<Variant> variants) {
        // TODO: this static method should not be filtering the unmapped reads.  To be addressed in another issue.
        JavaRDD<GATKRead> mappedReads = reads.filter(read -> ReadFilterLibrary.MAPPED.test(read));
        JavaPairRDD<GATKRead, Tuple2<Iterable<Variant>, ReferenceBases>> withVariantsWithRef;

        // Join Reads and Variants
        JavaPairRDD<GATKRead, Iterable<Variant>> withVariants = BroadcastJoinReadsWithVariants.join(mappedReads, variants);

        // Join Reads with ReferenceBases
        withVariantsWithRef = BroadcastJoinReadsWithRefBases.addBases(referenceDataflowSource, withVariants);
        return withVariantsWithRef.mapToPair(in -> new Tuple2<>(in._1(), new ReadContextData(in._2()._2(), in._2()._1())));
    }
}
