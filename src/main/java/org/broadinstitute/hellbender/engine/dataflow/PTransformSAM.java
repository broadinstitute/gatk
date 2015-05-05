package org.broadinstitute.hellbender.engine.dataflow;

import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * An abstract class representing a {@link PTransform<GATKRead,Output>}.
 * This provides a mechanism for propagating a {@link SAMFileHeader} to {@link DataFlowReadFn}..
 *
 * Calls to {@link #getHeader()} must be proceded by a call to {@link #setHeader(SAMFileHeader)}}
 * @param <Output> the output type of the resulting PCollection
 */
public abstract class PTransformSAM<Output> extends PTransform<PCollection<GATKRead>,PCollection<Output>> {
    private static final long serialVersionUID = 1l;

    private SAMFileHeader header;

    /**
     * @return a {@link SAMFileHeader}
     * @throws GATKException if {@link #setHeader(SAMFileHeader)} wasn't called before this.
     */
    public SAMFileHeader getHeader(){

        if (header == null){
            throw new GATKException("You must call setHeader before calling getHeader");
        }
        return header;
    }

    /**
     * @param header set the {@link SAMFileHeader} to provide to the {@link DataFlowReadFn} making up this transform.
     *                     Must be non-null and must be a valid header for the the PCollection<GATKRead> that this transform is processing.
     */
    public void setHeader(final SAMFileHeader header) {
        if (header == null){
            throw new IllegalArgumentException("null header");
        }
        this.header = header;
    }
}
