package org.broadinstitute.hellbender.engine.genomicsdb;


import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.*;
import org.broadinstitute.hellbender.cmdline.programgroups.GenomicsDBProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;


@CommandLineProgramProperties(
  summary = "Import streams of gVCFs GenomicsDB for given intervals",
  oneLineSummary = "Import gVCFs to GenomicsDB",
  programGroup = GenomicsDBProgramGroup.class
)

public class GenomicsDBWriter extends GATKTool {

  @ArgumentCollection
  public final VariantInputArgumentCollection variantArguments = requiresVariants() ?
    new RequiredVariantInputArgumentCollection() : new OptionalVariantInputArgumentCollection();

  @ArgumentCollection
  protected IntervalArgumentCollection intervalArgumentCollection = requiresIntervals() ? new RequiredIntervalArgumentCollection() : new OptionalIntervalArgumentCollection();

  /**
   * A complete traversal from start to finish. Tool authors who wish to "roll their own" traversal
   * from scratch can extend this class directly and implement this method. Walker authors should
   * instead extend a Walker class and implement the Walker-appropriate apply() method, since the
   * Walker base classes implement the various kinds of traversals for you.
   */
  @Override
  public void traverse() {

  }

  public boolean requiresVariants() {
    return true;
  }

}
