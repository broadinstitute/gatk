package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;

public interface MutableRead extends Read {
    void setName( String name );

    void setBases( byte[] bases );

    void setBaseQualities( byte[] baseQualities );

    void setPosition( String contig, int start );

    void setPosition( Locatable locatable );

    void setMatePosition( String contig, int start );

    void setMatePosition( Locatable locatable );

    void setFragmentLength( int fragmentLength );

    void setMappingQuality( int mappingQuality );

    void setCigar( Cigar cigar );

    void setCigar( String cigarString );

    void setReadGroup( String readGroupID );

    void setNumberOfReadsInFragment( int numberOfReads );

    void setReadNumber( int readNumber);

    void setIsPaired( boolean isPaired );

    void setIsProperPaired( boolean isProperPaired );

    void setIsUnmapped();

    void setMateIsUnmapped();

    void setIsNegativeStrand( boolean isNegativeStrand );

    void setMateIsNegativeStrand( boolean mateIsNegativeStrand );

    void setIsFirstOfPair( boolean isFirstOfPair );

    void setIsSecondOfPair( boolean isSecondOfPair );

    void setIsNonPrimaryAlignment( boolean isNonPrimaryAlignment );

    void setIsSupplementaryAlignment( boolean isSupplementaryAlignment );

    void setFailsVendorQualityCheck( boolean failsVendorQualityCheck );

    void setIsDuplicate( boolean isDuplicate );

    void setAttribute( String attributeName, Integer attributeValue );

    void setAttribute( String attributeName, String attributeValue );

    void setAttribute( String attributeName, byte[] attributeValue );

    void clearAttribute( String attributeName );

    void clearAttributes();
}