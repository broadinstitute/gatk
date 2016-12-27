package org.broadinstitute.hellbender.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class GenomicsDBProgramGroup implements CommandLineProgramGroup {
  @Override
  public String getName() { return "GenomicsDB"; }
  @Override
  public String getDescription() { return "Tools for writing and reading variants to and from GenomicsDB"; }
}