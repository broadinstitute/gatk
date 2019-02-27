package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Filter out reads that:
 *
 * <ul>
 *     <li>Fail platform/vendor quality checks (0x200)</li>
 *     <li>Are unmapped (0x4)</li>
 *     <li>Represent secondary/supplementary alignments (0x100 or 0x800)</li>
 * </ul>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY, summary = "Filter out reads that fail platform quality checks, are unmapped and represent secondary/supplementary alignments")
public final class MetricsReadFilter extends ReadFilter {
  // TODO: Should this be something more unique, such as a timestamp, in order
  // to behave with Spark's serialization?
  static final long serialVersionUID = 1L;

  private final boolean pfReadOnly;
  private final boolean alignedReadsOnly;

  // Command line parser requires a no-arg constructor
  public MetricsReadFilter() {this(true, true);}

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