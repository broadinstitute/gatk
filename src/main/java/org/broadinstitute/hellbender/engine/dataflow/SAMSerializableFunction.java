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

    private transient SAMFileHeader header;  //transient cache of the SAMFileHeader which is not Serializeable
    private final String headerString;
    private final SerializableFunction<SAMRecord, Output> f;

    /**
     * Wrap a function from SAMRecord -> Output into a function from Read -> Output
     * @param headerString A String to construct an appropriate SAMFileHeader from
     * @param f the SerializableFunction to wrap
     */
    public SAMSerializableFunction(String headerString, SerializableFunction<SAMRecord, Output> f){
        this.headerString = headerString;
        this.f = f;
    }

    private SAMFileHeader getHeader(){
        if (header == null) {
            final SAMTextHeaderCodec headerCodec = new SAMTextHeaderCodec();
            headerCodec.setValidationStringency(ValidationStringency.LENIENT);
            this.header = headerCodec.decode(new StringLineReader(headerString), "magic string");
        }
        return header;
    }


    @Override
    final public Output apply(Read read) {
        return f.apply(GenomicsConverter.makeSAMRecord(read, getHeader()));
    }
}
