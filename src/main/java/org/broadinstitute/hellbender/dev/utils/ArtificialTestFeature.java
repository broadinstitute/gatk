package org.broadinstitute.hellbender.dev.utils;

import htsjdk.tribble.Feature;

/**
 * A feature that's created programmatically.
 */
@SuppressWarnings("overrides")  // because I don't want to implement hashCode() but do need an equals() here
public class ArtificialTestFeature implements Feature {
  private String chr;
  private int start;
  private int end;

  public ArtificialTestFeature(final String chr, final int start, final int end) {
    this.chr = chr;
    this.start = start;
    this.end = end;
  }

  @Override
  public String getChr() {
    return chr;
  }

  @Override
  public String getContig() {
    return chr;
  }


  @Override
  public int getStart() {
    return start;
  }

  @Override
  public int getEnd() {
    return end;
  }

  @Override
  public boolean equals(Object other) {
    if (other == null || !(other instanceof ArtificialTestFeature)) {
      return false;
    }

    ArtificialTestFeature otherFeature = (ArtificialTestFeature) other;
    return chr.equals(otherFeature.getChr()) && start == otherFeature.getStart() && end == otherFeature.getEnd();
  }

  @Override
  public String toString() {
    return chr + ":" + start + "-" + end;   // (to improve output on test failures involving this class)
  }
}
