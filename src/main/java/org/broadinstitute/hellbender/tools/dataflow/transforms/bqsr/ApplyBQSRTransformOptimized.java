package org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.dataflow.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ContextShard;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsShard;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

/**
 * This transforms applies BQSR to the input reads, taking shards as input.
 */
public final class ApplyBQSRTransformOptimized extends PTransform<PCollection<ReadsShard>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1L;

    private final PCollectionView<SAMFileHeader> headerView;
    private final PCollectionView<BaseRecalOutput> recalView;
    private final ApplyBQSRArgumentCollection args;

    /**
     * @param headerView        The SAM header that corresponds to the reads you're going to pass as input.
     * @param recalibrationView the output from BaseRecalibration, with a single BaseRecalOutput object.
     * @param args              the recalibration args
     */
    public ApplyBQSRTransformOptimized(PCollectionView<SAMFileHeader> headerView, PCollectionView<BaseRecalOutput> recalibrationView, ApplyBQSRArgumentCollection args) {
        this.headerView = headerView;
        this.recalView = recalibrationView;
        this.args = args;
    }

    /**
     * @return The same reads as the input, but with updated quality scores.
     */
    @Override
    public PCollection<GATKRead> apply(PCollection<ReadsShard> input) {
        return input.apply(ParDo
            .named("ApplyBQSR")
            .withSideInputs(headerView, recalView)
            .of(new ApplyBQSR()));
    }

    private class ApplyBQSR extends DoFnWLog<ReadsShard, GATKRead> {
        private static final long serialVersionUID = 1L;
        private transient BQSRReadTransformer transformer;
        private transient SAMFileHeader header;

        public ApplyBQSR() {
            super("ApplyBQSR");
        }

        @Override
        public void processElement(ProcessContext c) {
            if (null==transformer) {
                header = c.sideInput(headerView);
                BaseRecalOutput info = c.sideInput(recalView);
                transformer = new BQSRReadTransformer(header, info, args.quantizationLevels, args.disableIndelQuals, args.PRESERVE_QSCORES_LESS_THAN, args.emitOriginalQuals, args.globalQScorePrior);
                // it's OK if this same object is used for multiple bundles
            }
            for (GATKRead r : c.element().reads) {
                setHeader(r, header);
                final GATKRead transformed = transformer.apply(r);
                setHeader(transformed, null);
                c.output(transformed);
            }
        }

        private void setHeader(GATKRead r, SAMFileHeader h) {
            if (r instanceof SAMRecordToGATKReadAdapter) {
                SAMRecordToGATKReadAdapter a = (SAMRecordToGATKReadAdapter)r;
                a.setHeader(h);
            }
        }

    }
}


