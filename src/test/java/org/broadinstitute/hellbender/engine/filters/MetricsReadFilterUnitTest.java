package org.broadinstitute.hellbender.engine.filters;

import static org.mockito.Mockito.*;
import static org.testng.Assert.*;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;

/**
 * No need to extend from GATKBaseTest here. Just testing some logic of {@link MetricsReadFilter}
 */
public class MetricsReadFilterUnitTest {
  private final MetricsReadFilter anyFilter =
      new MetricsReadFilter(false, false);
  private final MetricsReadFilter pfOnlyFilter =
      new MetricsReadFilter(true, false);
  private final MetricsReadFilter alignedOnlyFilter =
      new MetricsReadFilter(false, true);
  private final MetricsReadFilter pfAndAlignedOnlyFilter =
      new MetricsReadFilter(true, true);

  @Test
  public void blocksSecondaryAlignmentReads() {
    final GATKRead secondaryAlignedRead =
        when(mock(GATKRead.class).isSecondaryAlignment())
          .thenReturn(true).getMock();
    this.alwaysBlocks(secondaryAlignedRead);
  }

  @Test
  public void blocksSupplementaryAlignmentReads() {
    final GATKRead supplementaryAlignedRead =
        when(mock(GATKRead.class).isSupplementaryAlignment())
          .thenReturn(true).getMock();
    this.alwaysBlocks(supplementaryAlignedRead);
  }

  @Test
  public void canRestrictToPfReadsOnly() {
    final GATKRead pfRead = when(mock(GATKRead.class).failsVendorQualityCheck())
        .thenReturn(false).getMock();
    final GATKRead nonPfRead = when(
        mock(GATKRead.class).failsVendorQualityCheck())
          .thenReturn(true).getMock();
    assertTrue(this.pfOnlyFilter.test(pfRead), this.pfOnlyFilter.toString());
    assertFalse(this.pfOnlyFilter.test(nonPfRead),
        this.pfOnlyFilter.toString());
  }
  @Test
  public void canRestrictToAlignedReadsOnly() {
    final GATKRead alignedRead = when(mock(GATKRead.class).isUnmapped())
        .thenReturn(false).getMock();
    final GATKRead unalignedRead = when(
        mock(GATKRead.class).isUnmapped())
          .thenReturn(true).getMock();
    assertTrue(this.alignedOnlyFilter.test(alignedRead),
        this.alignedOnlyFilter.toString());
    assertFalse(this.alignedOnlyFilter.test(unalignedRead),
        this.alignedOnlyFilter.toString());
  }

  @Test
  public void canRestrictToBothAlignedAndPfReads() {
    final GATKRead alignedAndPfRead = mock(GATKRead.class);
    when(alignedAndPfRead.isUnmapped()).thenReturn(false);
    when(alignedAndPfRead.failsVendorQualityCheck()).thenReturn(false);
    assertTrue(this.pfAndAlignedOnlyFilter.test(alignedAndPfRead));

    final GATKRead alignedButNotPfRead = mock(GATKRead.class);
    when(alignedButNotPfRead.isUnmapped()).thenReturn(false);
    when(alignedButNotPfRead.failsVendorQualityCheck()).thenReturn(true);
    assertFalse(this.pfAndAlignedOnlyFilter.test(alignedButNotPfRead));

    final GATKRead pfButUnalignedRead = mock(GATKRead.class);
    when(pfButUnalignedRead.failsVendorQualityCheck()).thenReturn(false);
    when(pfButUnalignedRead.isUnmapped()).thenReturn(true);
    assertFalse(this.pfAndAlignedOnlyFilter.test(pfButUnalignedRead));
  }

  private Collection<MetricsReadFilter> allPossibleMetricsFilters() {
    final Collection<MetricsReadFilter> metricFilters = new ArrayList<>();
    metricFilters.add(this.anyFilter);
    metricFilters.add(this.pfOnlyFilter);
    metricFilters.add(this.alignedOnlyFilter);
    metricFilters.add(this.pfAndAlignedOnlyFilter);
    return metricFilters;
  }

  private void alwaysBlocks(final GATKRead read) {
    this.allPossibleMetricsFilters().forEach(f -> {
        assertFalse(f.test(read), f.toString());
    });
  }
}