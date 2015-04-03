package org.broadinstitute.hellbender.tools.dataflow;

import com.google.api.services.genomics.model.Read;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.Sum;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.readers.bam.ReadBAMTransform;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.Contig;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowSAMFn;
import org.broadinstitute.hellbender.engine.dataflow.DataflowTool;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(usage = "Count bases in dataflow", usageShort = "count bases", programGroup = DataFlowProgramGroup.class)
public class CountBasesDataflow extends DataflowTool{

  @Argument
  private String outputFile;

  @Argument
  private File clientSecret = new File("client_secret.json");

  @Argument
  String bam;

  private GenomicsFactory.OfflineAuth getAuth(GCSOptions options){
    try {
      return GCSOptions.Methods.createGCSAuth(options);
    } catch (IOException e) {
      throw new GATKException("Couldn't create a dataflow auth object.", e);
    }
  }

  private String getHeaderString(GCSOptions options, String bamPath) {
    try {
      Storage.Objects storageClient = GCSOptions.Methods.createStorageClient(options, getAuth(options));
      final BAMIO.ReaderAndIndex r = BAMIO.openBAMAndExposeIndex(storageClient, bamPath);
      SamReader reader = r.reader;
      return reader.getFileHeader().getTextHeader();
    } catch (IOException e) {
      throw new GATKException("Failed to read bam header from "+ bamPath+".", e);
    }
  }


  @Override
  protected void setupPipeline(Pipeline pipeline) {
    GCSOptions ops = PipelineOptionsFactory.fromArgs(new String[]{"--genomicsSecretsFile=" + clientSecret.getAbsolutePath()}).as(GCSOptions.class);
    GenomicsOptions.Methods.validateOptions(ops);


    GenomicsFactory.OfflineAuth auth = getAuth(ops);
    String headerString = getHeaderString(ops, bam);


    Iterable<Contig> contigs = intervals.stream()
            .map(i -> new Contig(i.getContig(), i.getStart(), i.getEnd()))
            .collect(Collectors.toList());


    PCollection<Read> preads= ReadBAMTransform.getReadsFromBAMFilesSharded(pipeline, auth, contigs, ImmutableList.of(bam));

    preads.apply(ParDo.of(new DataFlowSAMFn<Long>(headerString) {
      @Override
      protected void apply(SAMRecord read) {
        Long bases = (long) read.getReadBases().length;
        output(bases);
      }
    }))
      .apply(Sum.longsGlobally())
      .apply(ParDo.of(new DoFn<Long, String>() {
        @Override
        public void processElement(DoFn<Long, String>.ProcessContext c) throws Exception {
          c.output(String.valueOf(c.element()));
        }
      }))
      .apply(TextIO.Write.to(outputFile));
  }
}
