package org.broadinstitute.hellbender.dev.tools.walkers.bench;

import com.google.api.services.storage.Storage;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.util.Transport;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import com.google.genomics.v1.Read;

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
import java.io.InputStream;

/**
 * Reads the whole input. You probably don't want the Direct runner because then it tries to load the whole
 * thing in memory and goes "boom".
 */
@CommandLineProgramProperties(
        usage = "Times the reading of the input BAM's bytes.",
        usageShort = "Times the reading of the input BAM.",
        programGroup = ReadProgramGroup.class
)
public final class BenchReadBytesDataflow extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    private final static Logger logger = LogManager.getLogger(BenchReadBytesDataflow.class);

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
            .apply(ParDo.of(new DoFnWLog<Integer, Integer>("readBytes") {
              @Override
              public void processElement(ProcessContext c) throws Exception {
                final int MB = 1024 * 1024;
                String f = filename;
                InputStream inputStream = BucketUtils.openFile(f, c.getPipelineOptions());
                byte[] buf = new byte[MB];
                long total = 0;
                long next = 10 * MB;
                while (total < 300 * MB) {
                  total += inputStream.read(buf);
                  if (total > next) {
                    bunny.stepEnd("Read " + (total / MB) + " MB");
                    next += 10 * MB;
                  }
                }
              }
            }));
    }

    @Override
    protected void afterPipeline(Pipeline pipeline) {
    }

}
