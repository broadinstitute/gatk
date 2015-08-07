package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.ImmutableList;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.dev.BunnyLog;
import org.broadinstitute.hellbender.dev.DoFnWLog;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.DataflowReadFilter;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.coders.GATKReadCoder;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;


/**
 * Bug hunting. Run this with an input in the cloud to reproduce the bug.
 */
@CommandLineProgramProperties(
    summary = "Bug repro attempt.",
    oneLineSummary = "Bug repro attempt.",
    programGroup = ReadProgramGroup.class
)
public class BadTypeRepro extends DataflowCommandLineProgram {
  private static final long serialVersionUID = 1L;

  private static final Logger logger = LogManager.getLogger(BadTypeRepro.class);

  // ------------------------------------------
  // Command-line options


  @Argument(shortName = "I", fullName = "input")
  protected String input = null;

  // ------------------------------------------

  // the inputs to BQSR
  private SAMFileHeader header;

  private final BunnyLog bunny = new BunnyLog(logger);

  @Override
  protected void setupPipeline(Pipeline pipeline) {
    try {
      bunny.start("BaseRecalibratorDataflow");

      PCollection<GATKRead> reads = ingestReadsAndGrabHeader(pipeline, input);

    } catch (UserException | GATKException rx) {
      throw rx;
    } catch (Exception x) {
      throw new GATKException("Unexpected: " + x.getMessage(), x);
    }
  }

  /**
   * reads local disks or GCS -> header, and PCollection
   */
  private PCollection<GATKRead> ingestReadsAndGrabHeader(final Pipeline pipeline, String filename) throws IOException {
    String beforePath = filename;

    // input reads
    if (BucketUtils.isRemoteStorageUrl(beforePath)) {
      // set up ingestion on the cloud
      // but read the header locally
      InputStream inputstream = BucketUtils.openFile(beforePath, pipeline.getOptions());
      SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(inputstream));
      header = reader.getFileHeader();
      reader.close();

      final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
      final List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
      return new ReadsDataflowSource(beforePath, pipeline).getReadPCollection(intervals, ValidationStringency.SILENT);
    } else {
      // ingestion from local file
      System.out.println("To repro the bug, you're going to need to read from a GCS file.");
      try (ReadsDataSource readsSource = new ReadsDataSource(new File(beforePath))) {
        header = readsSource.getHeader();
        List<GATKRead> readLst = new ArrayList<>();
        ReadFilter readFilter = BaseRecalibratorWorker.readFilter(header);
        for (GATKRead read : readsSource) {
          if (readFilter.test(read)) {
            readLst.add(read);
          }
        }
        return pipeline.apply("input ingest", Create.of(readLst).withCoder(new GATKReadCoder()));
      }
    }
  }


}
