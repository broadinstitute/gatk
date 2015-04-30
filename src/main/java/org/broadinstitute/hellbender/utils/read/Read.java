package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;

public interface Read extends Locatable {
    String getName();

    int getLength();

    byte[] getBases();

    byte[] getBaseQualities();

    int getUnclippedStart();

    int getUnclippedEnd();

    String getMateContig();

    int getMateStart();

    int getFragmentLength();

    int getMappingQuality();

    Cigar getCigar();

    String getReadGroup();

    int getNumberOfReadsInFragment();

    int getReadNumber();

    boolean isPaired();

    boolean isProperlyPaired();

    boolean isUnmapped();

    boolean mateIsUnmapped();

    boolean isReverseStrand();

    boolean mateIsReverseStrand();

    boolean isFirstOfPair();

    boolean isSecondOfPair();

    boolean isNonPrimaryAlignment();

    boolean isSupplementaryAlignment();

    boolean failsVendorQualityCheck();

    boolean isDuplicate();

    boolean hasAttribute(String attributeName);

    Integer getAttributeAsInteger(String attributeName);

    String getAttributeAsString(String attributeName);

    byte[] getAttributeAsByteArray(String attributeName);

    Read copy();
}
