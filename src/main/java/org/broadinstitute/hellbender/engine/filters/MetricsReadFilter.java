package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.hellbender.utils.read.GATKRead;

public final class MetricsReadFilter implements ReadFilter {
  // TODO: Should this be something more unique, such as a timestamp, in order
  // to behave with Spark's serialization?
  static final long serialVersionUID = 1L;

  private final boolean pfReadOnly;
  private final boolean alignedReadsOnly;

  /**
   * @param alignedReadsOnly If set to true calculate mean quality over aligned
   * reads only.
   * @param pfReadOnly If set to true calculate mean quality over passing filter
   * reads only.
   */
  public MetricsReadFilter(final boolean pfReadOnly,
      final boolean alignedReadsOnly) {
    this.pfReadOnly = pfReadOnly;
    this.alignedReadsOnly = alignedReadsOnly;
  }
  @Override
  public boolean test(final GATKRead read) {
    return (!this.pfReadOnly || !read.failsVendorQualityCheck()) &&
        (!this.alignedReadsOnly || !read.isUnmapped()) &&
        !read.isSecondaryAlignment() && !read.isSupplementaryAlignment();
  }

  @Override
  public String toString() {
    return "[MetricsReadFilter{"
        + "pfReadOnly: " + this.pfReadOnly
        + ", alignedReadsOnly: " + this.alignedReadsOnly
        + "}]";
  }
}