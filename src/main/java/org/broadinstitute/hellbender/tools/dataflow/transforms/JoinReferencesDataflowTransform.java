package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.util.Map;

public class JoinReferencesDataflowTransform extends PTransformSAM<Long> {
    private static final long serialVersionUID = 1L;
    private static final Log log = Log.getInstance(JoinReferencesDataflowTransform.class);

    private Pipeline pipeline;
    private ReferenceFileSource refSource;

    public JoinReferencesDataflowTransform(Pipeline pipeline, ReferenceFileSource refSource) {
        this.pipeline = pipeline;
        this.refSource = refSource;
    }

    @Override
    public PCollection<Long> apply(final PCollection<GATKRead> reads) {
        PCollection<KV<String, ReferenceBases>> refCollection = null;
        try {
            log.info("Reading reference...");
            refCollection = pipeline.apply(Create.of(refSource.getAllReferenceBases()));
            log.info("... done reading reference");
        } catch (IOException e) {
            throw new IllegalStateException(e);
        }
        log.info("Creating refView...");
        final PCollectionView<Map<String, ReferenceBases>> refView = refCollection.apply(View.asMap());
        log.info("... done");
        return reads.apply(ParDo.of(new AddBasesFn(refView)).withSideInputs(refView))
                .apply(Count.globally());
    }

    private static class AddBasesFn extends DoFn<GATKRead, KV<GATKRead, ReferenceBases>> {
        private static final long serialVersionUID = 1L;

        private final PCollectionView<Map<String, ReferenceBases>> refView;

        public AddBasesFn(PCollectionView<Map<String, ReferenceBases>> refView) {
            this.refView = refView;
        }

        @Override
        public void processElement(ProcessContext c) throws Exception {
            GATKRead read = c.element();
            Map<String, ReferenceBases> ref = c.sideInput(refView);
            String refName = read.getContig();
            if (refName.startsWith("chr")) {
                refName = refName.substring(3);
            }
            int start = read.getStart();
            int end = read.getEnd();
            log.info("Read refName: " + refName);
            log.info("Read start: " + start);
            log.info("Read end: " + end);
            log.info("Ref names: " + ref.keySet());
            if (end < start) {
                log.info("End is less than start. Skipping.");
                return;
            }
            ReferenceBases bases = ref.get(refName);
            ReferenceBases subset = bases.getSubset(new SimpleInterval(refName, start, end));
            c.output(KV.of(read, subset));
        }
    }
}