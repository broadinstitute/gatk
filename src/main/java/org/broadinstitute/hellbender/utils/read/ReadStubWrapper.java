package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;

/**
 * Created by droazen on 5/15/15.
 */
public class ReadStubWrapper implements Read {

    private com.google.api.services.genomics.model.Read googleRead;

    public ReadStubWrapper( com.google.api.services.genomics.model.Read googleRead ) {
        this.googleRead = googleRead;
    }

    @Override
    public String getName() {
        return null;
    }

    @Override
    public int getLength() {
        return 0;
    }

    @Override
    public byte[] getBases() {
        return new byte[0];
    }

    @Override
    public byte[] getBaseQualities() {
        return new byte[0];
    }

    @Override
    public int getUnclippedStart() {
        return 0;
    }

    @Override
    public int getUnclippedEnd() {
        return 0;
    }

    @Override
    public String getMateContig() {
        return null;
    }

    @Override
    public int getMateStart() {
        return 0;
    }

    @Override
    public int getFragmentLength() {
        return 0;
    }

    @Override
    public int getMappingQuality() {
        return 0;
    }

    @Override
    public Cigar getCigar() {
        return null;
    }

    @Override
    public String getReadGroup() {
        return null;
    }

    @Override
    public int getNumberOfReadsInFragment() {
        return 0;
    }

    @Override
    public int getReadNumber() {
        return 0;
    }

    @Override
    public boolean isPaired() {
        return false;
    }

    @Override
    public boolean isProperlyPaired() {
        return false;
    }

    @Override
    public boolean isUnmapped() {
        return false;
    }

    @Override
    public boolean mateIsUnmapped() {
        return false;
    }

    @Override
    public boolean isReverseStrand() {
        return false;
    }

    @Override
    public boolean mateIsReverseStrand() {
        return false;
    }

    @Override
    public boolean isFirstOfPair() {
        return false;
    }

    @Override
    public boolean isSecondOfPair() {
        return false;
    }

    @Override
    public boolean isNonPrimaryAlignment() {
        return false;
    }

    @Override
    public boolean isSupplementaryAlignment() {
        return false;
    }

    @Override
    public boolean failsVendorQualityCheck() {
        return false;
    }

    @Override
    public boolean isDuplicate() {
        return false;
    }

    @Override
    public boolean hasAttribute( String attributeName ) {
        return false;
    }

    @Override
    public Integer getAttributeAsInteger( String attributeName ) {
        return null;
    }

    @Override
    public String getAttributeAsString( String attributeName ) {
        return null;
    }

    @Override
    public byte[] getAttributeAsByteArray( String attributeName ) {
        return new byte[0];
    }

    @Override
    public Read copy() {
        return null;
    }

    @Override
    public String getContig() {
        return null;
    }

    @Override
    public int getStart() {
        return 0;
    }

    @Override
    public int getEnd() {
        return 0;
    }
}
