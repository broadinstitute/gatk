package org.broadinstitute.hellbender.dev.tools.walkers.bench;

import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.ListBasesResponse;
import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.options.DataflowPipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.options.PipelineOptionsFactory;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.util.Transport;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.cloud.genomics.dataflow.utils.GenomicsOptions;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.base.Stopwatch;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.dev.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;

import java.io.IOException;
import java.security.GeneralSecurityException;
import java.util.concurrent.TimeUnit;

public final class BenchRefAPI {

  public static void main(String[] args) throws IOException {

    // hg19 is "referenceset EMWV_ZfLxrDY-wE"
    // that in turn contains stuff like reference "EMvD8Nr1mPeyAw"
    // GRCh38 is "referenceset EMud_c37lKPXTQ"
    // and that contains tons of stuff.
    // ELT13vHD8KjUEQ is chromosome 1 of GRCh38
    String refname = "ELT13vHD8KjUEQ";
    // where to start (keep in mind the first few are all 'N')
    long start = 1000000L;
    // how many bases to get.
    long count = 2864400L;

    DataflowPipelineOptions popts = PipelineOptionsFactory.fromArgs(args).create().as(DataflowPipelineOptions.class);
    popts.setSecretsFile("../client-secrets.json");
    Genomics genomicsService = createGenomicsService(popts);


    final Genomics.References.Bases.List listRequest = genomicsService.references().bases()
        .list(refname);
    listRequest.setStart(start);
    listRequest.setEnd(start+count);

    Stopwatch sw = Stopwatch.createStarted();
    ListBasesResponse bases = listRequest.execute();
    String acgtetc = bases.getSequence();
    sw.stop();

    System.out.println("Getting "+count+" bases from "+refname+" took "+sw.elapsed(TimeUnit.MILLISECONDS)+" ms.");

  }

  private static Genomics createGenomicsService( final PipelineOptions pipelineOptions ) {
    try {
      final GenomicsFactory.OfflineAuth auth = GenomicsOptions.Methods.getGenomicsAuth(pipelineOptions.as(GCSOptions.class));
      return auth.getGenomics(auth.getDefaultFactory());
    }
    catch ( GeneralSecurityException e ) {
      throw new UserException("Authentication failed for Google genomics service", e);
    }
    catch ( IOException e ) {
      throw new UserException("Unable to access Google genomics service", e);
    }
  }

}