package org.broadinstitute.hellbender.engine.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringLineReader;

import java.io.Serializable;
import java.util.function.Function;


public class SAMSerializableFunction<O> implements SerializableFunction<Read, O> {

        private transient SAMFileHeader header;
        private final String headerString;
        private final SerializableFunction<SAMRecord, O> f;

        public SAMSerializableFunction(String headerString, SerializableFunction<SAMRecord,O> f){
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
        final public O apply(Read read) {
            return f.apply(GenomicsConverter.makeSAMRecord(read, getHeader()));
        }
}
