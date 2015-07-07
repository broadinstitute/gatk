package org.broadinstitute.hellbender.engine.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;

/**
 * Transform that wraps Google Genomics {@link Read} objects within the {@link GATKRead} interface.
 *
 * No copy operations are performed -- this is simple encapsulation.
 */
public final class GoogleGenomicsReadToGATKRead extends PTransform<PCollection<Read>, PCollection<GATKRead>> {
    private static final long serialVersionUID = 1l;

    @Override
    public PCollection<GATKRead> apply( PCollection<Read> input ) {
        return input.apply(ParDo.of(new DoFn<Read, GATKRead>() {
            private static final long serialVersionUID = 1l;

            @Override
            public void processElement( ProcessContext c ) throws Exception {
                c.output(new GoogleGenomicsReadToGATKReadAdapter(c.element()));
            }
        }));
    }
}