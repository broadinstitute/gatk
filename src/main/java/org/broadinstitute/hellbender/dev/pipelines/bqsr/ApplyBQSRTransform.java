package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.dev.DoFnWLog;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalOutput;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * This transforms applies BQSR to the input reads.
 */
public final class ApplyBQSRTransform extends PTransform<PCollection<GATKRead>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1L;

    private final PCollectionView<SAMFileHeader> headerView;
    private final PCollectionView<BaseRecalOutput> recalView;
    private final ApplyBQSRArgumentCollection args;

    /**
     * @param headerView The SAM header that corresponds to the reads you're going to pass as input.
     * @param recalibrationView the output from BaseRecalibration, with a single BaseRecalOutput object.
     * @param args the recalibration args
     */
    public ApplyBQSRTransform(PCollectionView<SAMFileHeader> headerView, PCollectionView<BaseRecalOutput> recalibrationView, ApplyBQSRArgumentCollection args) {
        this.headerView = headerView;
        this.recalView = recalibrationView;
        this.args = args;
    }

    /**
     * @return The same reads as the input, but with updated quality scores.
     */
    @Override
    public PCollection<GATKRead> apply(PCollection<GATKRead> input) {
        return input.apply(ParDo
            .named("ApplyBQSR")
            .withSideInputs(headerView, recalView)
            .of(new ApplyBQSR()));
    }

    private class ApplyBQSR extends DoFnWLog<GATKRead, GATKRead> {
        private static final long serialVersionUID = 1L;
        private transient BQSRReadTransformer transformer;

        public ApplyBQSR() {
            super("ApplyBQSR");
        }

        @Override
        public void processElement(ProcessContext c) {
            if (null==transformer) {
                // set up the transformer, as needed
                SAMFileHeader header = c.sideInput(headerView);
                BaseRecalOutput info = c.sideInput(recalView);
                // it's OK if this same object is used for multiple bundles
                transformer = new BQSRReadTransformer(header, info, args.quantizationLevels, args.disableIndelQuals, args.PRESERVE_QSCORES_LESS_THAN, args.emitOriginalQuals, args.globalQScorePrior);
            }
            GATKRead r = c.element();
            final GATKRead transformed = transformer.apply(r);
            c.output(transformed);
        }

    }
}


