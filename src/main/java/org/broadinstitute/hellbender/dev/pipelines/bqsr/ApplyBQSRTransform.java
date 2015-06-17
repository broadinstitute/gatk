package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;

/**
 * This transforms applies BQSR to the input reads.
 */
public final class ApplyBQSRTransform extends PTransform<PCollection<Read>, PCollection<Read>> {
    private static final long serialVersionUID = 1L;

    private final SAMFileHeader header;
    private final PCollectionView<BaseRecalOutput> recalView;
    private final ApplyBQSRArgumentCollection args;

    /**
     * @param header The SAM header that corresponds to the reads you're going to pass as input.
     * @param recalibrationOutput the output from BaseRecalibration, with a single BaseRecalOutput object.
     * @param args the recalibration args
     */
    public ApplyBQSRTransform(SAMFileHeader header, PCollection<BaseRecalOutput> recalibrationOutput, ApplyBQSRArgumentCollection args) {
        this.header = header;
        this.recalView = recalibrationOutput.apply(View.asSingleton());
        this.args = args;
    }

    /**
     * @return The same reads as the input, but with updated quality scores.
     */
    @Override
    public PCollection<Read> apply(PCollection<Read> input) {
        return input.apply(ParDo
                .named("ApplyBQSR")
                .withSideInputs(recalView)
                .of(new ApplyBQSR()));
    }

    private class ApplyBQSR extends DoFn<Read, Read> {
        private static final long serialVersionUID = 1L;
        private final transient Logger logger = LogManager.getLogger(ApplyBQSR.class);
        private transient BQSRReadTransformer transformer;

        @Override
        public void processElement(ProcessContext c) {
            if (null==transformer) {
                // set up the transformer, as needed
                BaseRecalOutput info = c.sideInput(recalView);
                transformer = new BQSRReadTransformer(info, args.quantizationLevels, args.disableIndelQuals, args.PRESERVE_QSCORES_LESS_THAN, args.emitOriginalQuals, args.globalQScorePrior);
                // it's OK if this same object is used for multiple bundles
            }
            Read r = c.element();
            final SAMRecord sr = GenomicsConverter.makeSAMRecord(r, header);
            final SAMRecord transformed = transformer.apply(sr);
            try {
                Read e = ReadConverter.makeRead(transformed);
                c.output(e);
            } catch (SAMException x) {
                logger.warn("Skipping read " + sr.getReadName() + " because we can't convert it ("+x.getMessage()+").");
            } catch (NullPointerException y) {
                logger.warn("Skipping read " + sr.getReadName() + " because we can't convert it. (null?)");
            }
        }

    }
}


