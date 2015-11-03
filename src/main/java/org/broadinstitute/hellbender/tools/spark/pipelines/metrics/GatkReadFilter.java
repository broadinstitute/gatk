package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import org.apache.spark.api.java.function.Function;
import org.broadinstitute.hellbender.utils.read.GATKRead;

class GatkReadFilter {
  /** 
   * @param aligned_reads_only If set to true calculate mean quality over
   * aligned reads only.
   * @param pf_read_only If set to true calculate mean quality over passing
   * filter reads only.
   */
  public static Function<GATKRead, Boolean> by(final boolean pf_read_only,
      final boolean aligned_reads_only) {
    return read -> (!pf_read_only || !read.failsVendorQualityCheck()) &&
        (!aligned_reads_only || !read.isUnmapped()) &&
        !read.isSecondaryAlignment() &&
        !read.isSupplementaryAlignment();
  }
}