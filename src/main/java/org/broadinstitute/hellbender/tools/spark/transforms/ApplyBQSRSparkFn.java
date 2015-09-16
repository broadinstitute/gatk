package org.broadinstitute.hellbender.tools.spark.transforms;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;


public class ApplyBQSRSparkFn {

    public static JavaRDD<GATKRead> apply( JavaRDD<GATKRead> reads, final Broadcast<RecalibrationReport> reportBroadcast, final SAMFileHeader readsHeader ) {
        return reads.map(read -> {
            RecalibrationReport report = reportBroadcast.getValue();
            BQSRReadTransformer transformer = new BQSRReadTransformer(readsHeader, report, new ApplyBQSRArgumentCollection());
            // TODO: remove protective copy?
            return transformer.apply(read.copy());
        });
    }
}
