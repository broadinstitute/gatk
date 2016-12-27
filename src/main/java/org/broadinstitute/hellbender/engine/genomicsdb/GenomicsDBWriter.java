package org.broadinstitute.hellbender.engine.genomicsdb;


import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.GenomicsDBProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;

@CommandLineProgramProperties(
  summary = "Write stream of variants to GenomicsDB.",
  oneLineSummary = "Write variants to GenomicsDB",
  programGroup = GenomicsDBProgramGroup.class
)

public class GenomicsDBWriter extends GATKTool {
  /**
   * A complete traversal from start to finish. Tool authors who wish to "roll their own" traversal
   * from scratch can extend this class directly and implement this method. Walker authors should
   * instead extend a Walker class and implement the Walker-appropriate apply() method, since the
   * Walker base classes implement the various kinds of traversals for you.
   */
  @Override
  public void traverse() {

  }

}
