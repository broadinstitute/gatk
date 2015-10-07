package org.broadinstitute.hellbender.engine.spark;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceDataflowSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
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
 * inside the {@link ReferenceDataflowSource} argument to {@link #add}. See
 * {@link org.broadinstitute.hellbender.engine.dataflow.datasources.RefWindowFunctions} for examples.
 */
public class AddContextDataToReadSpark {
    public static JavaPairRDD<GATKRead, ReadContextData> add(
            final JavaRDD<GATKRead> reads, final ReferenceDataflowSource referenceDataflowSource,
            final JavaRDD<Variant> variants, final JoinStrategy joinStrategy) {
        // TODO: this static method should not be filtering the unmapped reads.  To be addressed in another issue.
        JavaRDD<GATKRead> mappedReads = reads.filter(read -> ReadFilterLibrary.MAPPED.test(read));
        // Join Reads and Variants
        JavaPairRDD<GATKRead, Iterable<Variant>> withVariants = JoinReadsWithVariants.join(mappedReads, variants);
        // Join Reads with ReferenceBases
        JavaPairRDD<GATKRead, Tuple2<Iterable<Variant>, ReferenceBases>> withVariantsWithRef;
        if (joinStrategy.equals(JoinStrategy.BROADCAST)) {
            withVariantsWithRef = BroadcastJoinReadsWithRefBases.addBases(referenceDataflowSource, withVariants);
        } else if (joinStrategy.equals(JoinStrategy.SHUFFLE)) {
            withVariantsWithRef = ShuffleJoinReadsWithRefBases.addBases(referenceDataflowSource, withVariants);
        } else {
            throw new UserException("Unknown JoinStrategy");
        }
        return withVariantsWithRef.mapToPair(in -> new Tuple2<>(in._1(), new ReadContextData(in._2()._2(), in._2()._1())));
    }
}
