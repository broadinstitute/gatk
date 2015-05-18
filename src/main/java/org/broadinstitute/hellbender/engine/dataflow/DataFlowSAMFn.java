package org.broadinstitute.hellbender.engine.dataflow;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

/**
 * A DoFn from Read -> Something which presents a SAMRecord to its apply method.
 * This is a temporary workaround for not having a common Read interface.
 */
public abstract class DataFlowSAMFn<Output> extends DoFn<Read, Output> {
    private static final long serialVersionUID = 1l;
    /**
     * field to cache the reconstructed SAMFileHeader
     * SAMFileHeader is not Serializable, but this class must be, so this is marked transient and reconstructed as needed
     */
    private final SAMFileHeader header;

    private final List<Output> toEmit = new ArrayList<>();

    /**
     * Construct a new DataFlowSAMFn and give it an appropriate {@link SAMFileHeader}
     * @param header A header String from a SAM or BAM file.  This can be read from a SAM or BAM file using {@link SAMFileHeader#getTextHeader()}
     */
    public DataFlowSAMFn(final SAMFileHeader header){
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
        final SAMRecord sam = GenomicsConverter.makeSAMRecord(c.element(), getHeader());
        apply(sam);
        toEmit.forEach(c::output);
        toEmit.clear();
    }

    /**
     * Subclasses must call this from within their {@link #apply(SAMRecord)} call, in order to output results.
     * May be called any number of times to output multiple values.
     * This should only be called from within apply, or the results will be undetermined.
     * @param output a value to output
     */
    protected final void output(Output output){
        toEmit.add(output);
    }

    /**
     * Subclasses must override apply.  It will be called once for every {@link Read} in the input PCollection.
     * Output from apply should be provide by calling {@link #output(Object)}
     * @param read a read to process
     */
    protected abstract void apply(SAMRecord read);
}
