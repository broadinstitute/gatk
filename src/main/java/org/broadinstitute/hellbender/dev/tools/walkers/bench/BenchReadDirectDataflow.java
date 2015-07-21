package org.broadinstitute.hellbender.dev.tools.walkers.bench;

import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.util.Transport;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;

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

/**
 * Reads the whole input. You probably don't want the Direct runner because then it tries to load the whole
 * thing in memory and goes "boom".
 */
@CommandLineProgramProperties(
        usage = "Times the reading of the input BAM and parsing to SAMRecord.",
        usageShort = "Times the parsing of the input BAM.",
        programGroup = ReadProgramGroup.class
)
public final class BenchReadDirectDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    private final static Logger logger = LogManager.getLogger(BenchReadDirectDataflow.class);

    @ArgumentCollection
    public final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected void setupPipeline(Pipeline pipeline) {
        if (readArguments.getReadFilesNames().size()>1) {
            throw new UserException("Sorry, we only support a single input file for now.");
        }
        final String filename = readArguments.getReadFilesNames().get(0);

      if (BucketUtils.fileExists(filename, pipeline.getOptions())) {
        System.out.println("File '"+filename+"' exists, OK");
      } else {
        throw new UserException.CouldNotReadInputFile(filename);
      }

      pipeline.apply(Create.of(1))
          .apply(ParDo.of(new DoFnWLog<Integer, Integer>("readDirect") {
            @Override
            public void processElement(ProcessContext c) throws Exception {
              final int MB = 1024 * 1024;
              String f = filename;

              Storage.Objects storage = Transport.newStorageClient(c.getPipelineOptions().as(GCSOptions.class)).build().objects();
              SamReader reader = BAMIO.openBAM(storage, filename, ValidationStringency.SILENT);

              SAMRecordIterator it = reader.iterator();
              int next = 1000000;
              // that should be about 300MB of reads
              for (int i=0; i<2864400; i++) {
                SAMRecord nextRecord = it.next();
                if (i==next) {
                  bunny.stepEnd("" + i + " reads");
                  next += 1000000;
                }
              }
            }}));
    }

    @Override
    protected void afterPipeline(Pipeline pipeline) {
    }

}
