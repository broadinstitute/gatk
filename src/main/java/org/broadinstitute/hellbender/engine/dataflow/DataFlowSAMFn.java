package org.broadinstitute.hellbender.engine.dataflow;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringLineReader;

import java.util.ArrayList;
import java.util.List;

public abstract class DataFlowSAMFn<O> extends DoFn<Read, O> {

  private transient SAMFileHeader header;
  private final String headerString;

  private List<O> toEmit = new ArrayList<>();

  public DataFlowSAMFn(String headerString){
    this.headerString = headerString;
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
  public void processElement(ProcessContext c) throws Exception {
    SAMRecord sam = GenomicsConverter.makeSAMRecord(c.element(), getHeader());
    apply(sam);
    toEmit.forEach(c::output);
    toEmit.clear();
  }


  protected void output(O output){
    toEmit.add(output);
  }

  protected abstract void apply(SAMRecord read);
}
