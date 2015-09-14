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
            ApplyBQSRArgumentCollection args = new ApplyBQSRArgumentCollection();
            BQSRReadTransformer transformer = new BQSRReadTransformer(readsHeader, report, args.quantizationLevels, args.disableIndelQuals, args.PRESERVE_QSCORES_LESS_THAN, args.emitOriginalQuals, args.globalQScorePrior);
            // TODO: remove protective copy?
            return transformer.apply(read.copy());
        });
    }
}
