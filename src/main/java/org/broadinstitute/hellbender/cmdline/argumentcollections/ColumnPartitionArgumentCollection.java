package org.broadinstitute.hellbender.cmdline.argumentcollections;


import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public final class ColumnPartitionArgumentCollection implements Serializable {
  private static final long serialVersionUID = 1L;

  @Argument(fullName = "columnPartitions",
    shortName = "CL",
    doc = "One or more genomic intervals over which to partition the GenomicsDB data",
    optional = true)
  protected final List<Long> columnPartitions = new ArrayList<>();

  protected List<Long> getPartitions() {
    return columnPartitions;
  }
  
}
