package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;

import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class ApplyBQSRSparkFn {

    public static JavaRDD<GATKRead> apply(JavaRDD<GATKRead> reads, final Broadcast<RecalibrationReport> reportBroadcast, final SAMFileHeader readsHeader, ApplyBQSRArgumentCollection args) {
        return reads.mapPartitions(readsIterator -> {
            final RecalibrationReport report = reportBroadcast.getValue();
            final BQSRReadTransformer transformer = new BQSRReadTransformer(readsHeader, report, args);//reuse this for all reads in the partition
            final Iterable<GATKRead> readsIterable = () -> readsIterator;
            return Utils.stream(readsIterable).map(read -> transformer.apply(read)).collect(Collectors.toList()).iterator();
        });
    }
}
