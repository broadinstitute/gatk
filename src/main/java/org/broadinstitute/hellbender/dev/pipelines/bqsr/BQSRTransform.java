package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionTuple;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;

/**
 *
 * BaseRecalibration, presented as a PTransform.
 *
 * Create the inputs like this:
 *
 * PCollection<Read> reads = ...
 * PCollection<SimpleInterval> knownIntervals = ...
 * PCollectionTuple inputs =
 *     PCollectionTuple.of(BaseRecalibratorDataflowUtils.readTag, reads)
 *                     .and(BaseRecalibratorDataflowUtils.intervalTag, intervals);
 */
public final class BQSRTransform extends PTransform<PCollectionTuple, PCollection<BaseRecalOutput>> {

    private final SAMFileHeader header;
    // local or GCS
    private final String referencePath;
    private final BaseRecalibrationArgumentCollection toolArgs;

    /**
     * @param header The SAM header that corresponds to the reads you're going to pass as input.
     * @param referencePath A path (local or GCS) to the reference file. It must also have a .dict and a .fasta.fai next to it.
     * @param toolArgs User-provided arguments that configure BaseRecalibration.
     */
    public BQSRTransform(SAMFileHeader header, String referencePath, BaseRecalibrationArgumentCollection toolArgs) {
        this.header = header;
        this.referencePath = referencePath;
        this.toolArgs = toolArgs;
    }

    /**
     * @return a single RecalibrationTables object that represents the output of phase 1 of BQSR.
     */
    @Override
    public PCollection<BaseRecalOutput> apply(PCollectionTuple input) {
        return BaseRecalibratorDataflowUtils.getRecalibrationOutput(header, input.get(BaseRecalibratorDataflowUtils.readTag), referencePath, toolArgs, input.get(BaseRecalibratorDataflowUtils.intervalTag));
    }


}
