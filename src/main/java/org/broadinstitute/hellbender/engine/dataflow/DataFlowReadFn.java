package org.broadinstitute.hellbender.engine.dataflow;


import com.google.cloud.dataflow.sdk.transforms.DoFn;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

/**
 * A DoFn from GATKRead -> Something which presents a GATKRead to its apply method.
 */
public abstract class DataFlowReadFn<Output> extends DoFn<GATKRead, Output> {
    private static final long serialVersionUID = 1l;
    /**
     * field to cache the SAMFileHeader
     */
    private final SAMFileHeader header;

    private final List<Output> toEmit = new ArrayList<>();

    /**
     * Construct a new DataFlowSAMFn and give it an appropriate {@link SAMFileHeader}
     * @param header A header String from a SAM or BAM file.  This can be read from a SAM or BAM file using {@link SAMFileHeader#getTextHeader()}
     */
    public DataFlowReadFn( final SAMFileHeader header ){
        this.header = header;
    }

    /**
     * @return a SAMFileHeader reconstituted from the header string this was constructed with
     */
    public final SAMFileHeader getHeader(){
        return header;
    }


    @Override
    public final void processElement(ProcessContext c) throws Exception {
        toEmit.clear();
        apply(c.element());
        toEmit.forEach(c::output);
        toEmit.clear();
    }

    /**
     * Subclasses must call this from within their {@link #apply(GATKRead)} call, in order to output results.
     * May be called any number of times to output multiple values.
     * This should only be called from within apply, or the results will be undetermined.
     * @param output a value to output
     */
    protected final void output(Output output){
        toEmit.add(output);
    }

    /**
     * Subclasses must override apply.  It will be called once for every {@link GATKRead} in the input PCollection.
     * Output from apply should be provide by calling {@link #output(Object)}
     * @param read a read to process
     */
    protected abstract void apply(GATKRead read);
}
