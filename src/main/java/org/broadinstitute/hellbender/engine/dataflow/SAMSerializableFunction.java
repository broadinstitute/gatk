package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringLineReader;


/**
 * A wrapper for converting {@link SerializableFunction<SAMRecord,Output>} into a {@link SerializableFunction<Read,Output>}
 * Provides automatic conversion from Read to SAMRecord
 * @param <Output>
 */
public final class SAMSerializableFunction<Output> implements SerializableFunction<Read, Output> {
    private static final long serialVersionUID = 1L;

    private SAMFileHeader header;
    private final SerializableFunction<SAMRecord, Output> f;

    /**
     * Wrap a function from SAMRecord -> Output into a function from Read -> Output
     * @param header A {@link SAMFileHeader } that is associated with the Reads this will be used on.
     * @param f the SerializableFunction to wrap
     */
    public SAMSerializableFunction(SAMFileHeader header, SerializableFunction<SAMRecord, Output> f){
        this.header = header;
        this.f = f;
    }

    private SAMFileHeader getHeader(){
        return header;
    }


    @Override
    final public Output apply(Read read) {
        return f.apply(GenomicsConverter.makeSAMRecord(read, getHeader()));
    }
}
