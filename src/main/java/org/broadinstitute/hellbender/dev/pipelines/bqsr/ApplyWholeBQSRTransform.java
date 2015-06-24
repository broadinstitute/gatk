package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionTuple;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;

/**
 * Both BQSR phases, together.
 *
 * Create the inputs like this:
 *
 * PCollection<Read> reads = ...
 * PCollection<SimpleInterval> knownIntervals = ...
 * PCollectionTuple inputs =
 *     PCollectionTuple.of(BaseRecalibratorDataflowUtils.readTag, reads)
 *                     .and(BaseRecalibratorDataflowUtils.intervalTag, intervals);
 */
public class ApplyWholeBQSRTransform extends PTransform<PCollectionTuple, PCollection<Read>> {


    private final SAMFileHeader header;
    // local or GCS
    private final String referencePath;
    private final BaseRecalibrationArgumentCollection recalArgs;
    private final ApplyBQSRArgumentCollection applyArgs;

    /**
     * @param header The SAM header that corresponds to the reads you're going to pass as input.
     * @param referencePath A path (local or GCS) to the reference file. It must also have a .dict and a .fasta.fai next to it.
     * @param recalArgs User-provided arguments that configure BaseRecalibration.
     * @param applyArgs User-provided arguments that configure ApplyBQSR.
     */
    public ApplyWholeBQSRTransform(SAMFileHeader header, String referencePath, BaseRecalibrationArgumentCollection recalArgs, ApplyBQSRArgumentCollection applyArgs) {
        this.header = header;
        this.referencePath = referencePath;
        this.recalArgs = recalArgs;
        this.applyArgs = applyArgs;
    }

    /**
     * @return The same reads as the input, but with updated quality scores.
     */
    @Override
    public PCollection<Read> apply(PCollectionTuple input) {
        PCollection<BaseRecalOutput> step1 = input.apply(new BQSRTransform(header, referencePath, recalArgs));
        ApplyBQSRTransform step2 = new ApplyBQSRTransform(header, step1, applyArgs);
        return input.get(BaseRecalibratorDataflowUtils.readTag).apply(step2);
    }
}
