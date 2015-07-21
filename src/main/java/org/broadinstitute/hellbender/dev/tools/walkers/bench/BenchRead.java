package org.broadinstitute.hellbender.dev.tools.walkers.bench;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.values.PCollection;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.dev.BunnyLog;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

@CommandLineProgramProperties(
        usage = "Times the loading and parsing of the input BAM.",
        usageShort = "Times the parsing of the input BAM.",
        programGroup = ReadProgramGroup.class
)
public final class BenchRead extends ReadWalker {
    private static final long serialVersionUID = 1L;
    private final static Logger logger = LogManager.getLogger(BenchRead.class);
  private final BunnyLog bunny = new BunnyLog(logger);
  private long count = 0;

  @Override
  public void onTraversalStart() {
    bunny.start("BenchRead");
  }

  @Override
  public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
    // yay I do nothing
    count++;
  }

  @Override
  public Object onTraversalDone() {
    logger.info("Processed "+count+" reads");
    bunny.end();
    return super.onTraversalDone();
  }

}
