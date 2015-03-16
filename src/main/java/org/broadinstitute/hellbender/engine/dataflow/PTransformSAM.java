package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.PTransform;
import com.google.cloud.dataflow.sdk.values.PCollection;
import org.broadinstitute.hellbender.exceptions.GATKException;

/**
 * An abstract class representing a {@link PTransform<Read, Output>}.
 * This provides a mechanism for propagating a {@link htsjdk.samtools.SAMFileHeader} to {@link DataFlowSAMFn}.
 * This is a workaround the fact that SAMFileHeader is not Serializable.
 *
 * Calls to {@link #getHeaderString()) must be proceded by a call to {@link #SetHeaderString())
 * @param <Output> the output type of the resulting PCollection
 */
public abstract class PTransformSAM<Output> extends PTransform<PCollection<Read>,PCollection<Output>> {
    private String headerString;

    /**
     * @return a String representation of a {@link htsjdk.samtools.SAMFileHeader}.
     * @throws GATKException if {@link #setHeaderString(String)} wasn't called before this.
     */
    public String getHeaderString(){
        if (headerString == null){
            throw new GATKException("You must call setHeaderString before calling getHeaderString");
        }
        return headerString;
    }

    /**
     * @param headerString set the headerString to provide to the {@link DataFlowSAMFn} making up this transform.
     *                     Must be non-null and must be a valid header for the the PCollection<Read> that this transform is processing.
     */
    public void setHeaderString(String headerString) {
        this.headerString = headerString;
    }
}
