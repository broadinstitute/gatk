package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * A simple filter for reads on dataflow
 */
public final class DataflowReadFilter extends PTransform<PCollection<GATKRead>, PCollection<GATKRead>> implements Serializable {
    private static final long serialVersionUID = 1L;
    
    private final ReadFilter readFilter;
    private final SAMFileHeader header;

    /**
     * ReadsFilter will use the given header to interpret each read, and then let the Read through
     * if the passed filter accepts it.
     */
    public DataflowReadFilter( final ReadFilter filterToApply, final SAMFileHeader header ) {
        Utils.nonNull(filterToApply);
        Utils.nonNull(header);

        this.readFilter = filterToApply;
        this.header = header;
    }

    /**
     * Filter out reads we don't want.
     */
    @Override
    public PCollection<GATKRead> apply(PCollection<GATKRead> in) {
        return in.apply(ParDo
                .named(readFilter.getClass().getSimpleName())
                .of(new DoFn<GATKRead, GATKRead>() {
                    private static final long serialVersionUID = 1L;
                    @Override
                    public void processElement(DoFn<GATKRead,GATKRead>.ProcessContext c) throws Exception {
                        GATKRead read = c.element();
                        if (readFilter.test(read)) {
                            c.output(read);
                        }
                    }
                }));
    }


}
