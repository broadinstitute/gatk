package org.broadinstitute.hellbender.utils.read;

import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.UUID;

/**
 * Unified read interface for use throughout the GATK.
 *
 * Adapter classes implementing this interface exist for both htsjdk's {@link SAMRecord} ({@link SAMRecordToGATKReadAdapter})
 * and the Google Genomics {@link Read} ({@link GoogleGenomicsReadToGATKReadAdapter})
 *
 * Since the adapter classes wrap the raw reads without making a copy, care must be taken to avoid
 * exposing the encapsulated reads, particularly if modifying reads in-place. As a result, this interface
 * should probably only be implemented by core engine-level classes.
 *
 * All GATKRead methods that return mutable reference types make defensive copies, with the exception
 * of the conversion methods {@link #convertToSAMRecord} and {@link #convertToGoogleGenomicsRead}.
 */
public interface GATKRead extends Locatable {

    /**
     * @return A globally unique identifier for this read, suitable for use as a key
     */
    UUID getUUID();

    /**
     * @return The name of the read (equivalent to QNAME in SAM), or {@code null} if the read has no name.
     */
    String getName();

    /**
     * Set the name of the read (equivalent to QNAME in SAM), or set to {@code null} if the read has no name.
     *
     * @param name new name for the read
     */
    void setName( final String name );

    /**
     * @return The number of bases in the read
     *
     * Note: This is not necessarily the same as the number of reference bases the read is aligned to.
     */
    default int getLength() {
        return getBases().length;
    }

    /**
     * @return True if the read has no bases, otherwise false
     */
    default boolean isEmpty() {
        return getLength() == 0;
    }

    /**
     * Set the position of the read (the contig and the start position). Cannot be used to
     * set the read to an unmapped position; use {@link #setIsUnmapped} for that purpose.
     *
     * @param contig Contig the read is mapped to
     * @param start Start position of the read (1-based, inclusive)
     * @throws IllegalArgumentException if contig is null or "*", or start is less than 1
     */
    void setPosition( final String contig, final int start );

    /**
     * Set the position of the read using the position of an existing {@link Locatable}. Cannot be used to
     * set the read to an unmapped position; use {@link #setIsUnmapped} for that purpose.
     *
     * @param locatable {@link Locatable} representing the 1-based, inclusive position to assign to the read
     * @throws IllegalArgumentException if locatable is null, or its contig is null or "*", or its start position is less than 1
     */
    void setPosition( final Locatable locatable );

    /**
     * Returns the alignment start (1-based, inclusive) adjusted for clipped bases.
     * For example, if the read has an alignment start of 100 but the first 4 bases
     * were clipped (hard or soft clipped) then this method will return 96.
     *
     * For unmapped reads, always returns {@link ReadConstants#UNSET_POSITION}
     *
     * @return The alignment start (1-based, inclusive) adjusted for clipped bases,
     *         or {@link ReadConstants#UNSET_POSITION} if the read is unmapped.
     */
    int getUnclippedStart();

    /**
     * Returns the alignment end (1-based, inclusive) adjusted for clipped bases.
     * For example, if the read has an alignment end of 100 but the last 7 bases
     * were clipped (hard or soft clipped) then this method will return 107.
     *
     * For unmapped reads, always returns {@link ReadConstants#UNSET_POSITION}
     *
     * @return The alignment end (1-based, inclusive) adjusted for clipped bases,
     *         or {@link ReadConstants#UNSET_POSITION} if the read is unmapped.
     */
    int getUnclippedEnd();

    /**
     * @return The contig that this read's mate is mapped to, or {@code null} if the mate is unmapped
     * @throws IllegalStateException if the read is not paired (has no mate)
     */
    String getMateContig();

    /**
     * @return The alignment start (1-based, inclusive) of this read's mate, or {@link ReadConstants#UNSET_POSITION}
     *         if the mate is unmapped.
     * @throws IllegalStateException if the read is not paired (has no mate)
     */
    int getMateStart();

    /**
     * Set the position of the read's mate (the contig and the start position). Cannot be used to
     * set the read's mate to an unmapped position; use {@link #setMateIsUnmapped} for that purpose.
     *
     * Calling this method has the additional effect of marking the read as paired, as if {@link #setIsPaired}
     * were invoked with true.
     *
     * @param contig Contig the read's mate is mapped to
     * @param start Start position of the read's mate (1-based, inclusive)
     * @throws IllegalArgumentException if contig is null or "*", or start is less than 1
     */
    void setMatePosition( final String contig, final int start );

    /**
     * Set the position of the read's mate using the position of an existing {@link Locatable}. Cannot be used to
     * set the read's mate to an unmapped position; use {@link #setMateIsUnmapped} for that purpose.
     *
     * Calling this method has the additional effect of marking the read as paired, as if {@link #setIsPaired}
     * were invoked with true.
     *
     * @param locatable {@link Locatable} representing the 1-based, inclusive position to assign to the read's mate
     * @throws IllegalArgumentException if locatable is null, or its contig is null or "*", or its start position is less than 1
     */
    void setMatePosition( final Locatable locatable );

    /**
     * Returns the observed length of the read's fragment (equivalent to TLEN in SAM).
     *
     * Warning: the precise meaning of this field is implementation/technology dependent.
     *
     * @return The observed length of the fragment (equivalent to TLEN in SAM), or 0 if unknown.
     *         Negative if the mate maps to a lower position than the read.
     */
    int getFragmentLength();

    /**
     * Set the observed length of the read's fragment (equivalent to TLEN in SAM).
     *
     * Warning: the precise meaning of this field is implementation/technology dependent.
     *
     * @param fragmentLength Observed length of the read's fragment; may be negative
     *                       if the mate maps to a lower position than the read,
     *                       or 0 if unknown.
     */
    void setFragmentLength( final int fragmentLength );

    /**
     * @return The mapping quality of this alignment, representing the phred-scaled likelihood that the read maps
     *         to this position as opposed to other locations. Returns {@link ReadConstants#NO_MAPPING_QUALITY}
     *         if there is none.
     */
    int getMappingQuality();

    /**
     * Set the mapping quality of this alignment, representing how likely the read maps to this position as
     * opposed to other locations. Set to {@link ReadConstants#NO_MAPPING_QUALITY} if there is none.
     *
     * @param mappingQuality mapping quality of this alignment; must be between 0 and 255, inclusive
     * @throws IllegalArgumentException if the mapping quality is less than 0 or greater than 255
     */
    void setMappingQuality( final int mappingQuality );

    /**
     * @return The read sequence as ASCII bytes ACGTN=, or an empty byte[] if no sequence is present.
     *
     * This method makes a defensive copy of the bases array before returning it, so modifying the
     * returned array will not alter the bases in the read.
     */
    byte[] getBases();

    /**
     * @return All bases in the read as a single String, or {@link ReadConstants#NULL_SEQUENCE_STRING}
     *         if the read is empty.
     */
    default String getBasesString() {
        return isEmpty() ? ReadConstants.NULL_SEQUENCE_STRING : StringUtil.bytesToString(getBases());
    }

    /**
     * Set the read's sequence.
     *
     * @param bases The read sequence as ASCII bytes ACGTN=. May be empty or null if no sequence is present.
     */
    void setBases( final byte[] bases );

    /**
     * @return Base qualities as binary phred scores (not ASCII), or an empty byte[] if base qualities are not present.
     *
     * This method makes a defensive copy of the base qualities array before returning it, so modifying the
     * returned array will not alter the base qualities in the read.
     */
    byte[] getBaseQualities();

    /**
     * Set the read's base qualities.
     *
     * @param baseQualities Base qualities as binary phred scores (not ASCII); negative values not allowed.
     *                      May be empty or null if no base qualities are present.
     * @throws GATKException if an invalid (negative) base quality is provided
     */
    void setBaseQualities( final byte[] baseQualities );

    /**
     * @return Cigar object describing how the read aligns to the reference, or an empty Cigar object if no cigar is present.
     *
     * This method makes a defensive copy of the Cigar within the read if necessary, so modifying the return value of
     * this method will not modify the read's Cigar.
     */
    Cigar getCigar();

    /**
     * Set the read's Cigar using an existing {@link Cigar} object describing how the read aligns to the reference.
     *
     * @param cigar {@link Cigar} object describing how the read aligns to the reference; May be null or empty
     *              if the read has none.
     */
    void setCigar( final Cigar cigar );

    /**
     * Set the read's Cigar using a textual cigar String describing how the read aligns to the reference.
     *
     * @param cigarString Cigar String describing how the read aligns to the reference; May be null or empty
     *                    if the read has none.
     */
    void setCigar( final String cigarString );

    /**
     * @return The ID of the read group this read belongs to, or {@code null} for none.
     */
    String getReadGroup();

    /**
     * Set the ID of the read group this read belongs to, or {@code null} for none.
     *
     * @param readGroupID ID of the read group this read belongs to, or {@code null} for none.
     */
    void setReadGroup( final String readGroupID );

    /**
     * @return True if this read is paired (ie., has a mate), otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isPaired();

    /**
     * Mark the read as paired (having a mate) or not paired.
     *
     * Setting this to false has the additional effect of marking the read as not
     * properly paired, as if {@link #setIsProperlyPaired} were invoked with false.
     *
     * @param isPaired True if this read is paired (ie., has a mate), otherwise false.
     */
    void setIsPaired( final boolean isPaired );

    /**
     * @return True if this read is paired and the orientation and the distance between reads from the fragment are
     *         consistent with the sequencing protocol, otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isProperlyPaired();

    /**
     * Mark the read as properly paired (or not properly paired).
     *
     * Setting this to true has the additional effect of marking the read as paired,
     * as if {@link #setIsPaired} were invoked with true.
     *
     * @param isProperlyPaired True if this read is paired and the orientation and the distance between reads from
     *                         the fragment are consistent with the sequencing protocol, otherwise false.
     */
    void setIsProperlyPaired( final boolean isProperlyPaired );

    /**
     * @return True if this read is unmapped (this includes reads that have a position but are explicitly marked as unmapped,
     *         as well as reads that lack a fully-defined position but are not explicitly marked as unmapped). Otherwise false.
     */
    boolean isUnmapped();

    /**
     * Mark the read as unmapped (lacking a defined position on the genome).
     *
     * To mark a read as mapped, use {@link #setPosition}
     */
    void setIsUnmapped();

    /**
     * @return True if this read's mate is unmapped (this includes mates that have a position but are explicitly marked as unmapped,
     *         as well as mates that lack a fully-defined position but are not explicitly marked as unmapped). Otherwise false.
     * @throws IllegalStateException if the read is not paired (has no mate)
     */
    boolean mateIsUnmapped();

    /**
     * Mark the read's mate as unmapped (lacking a defined position on the genome).
     *
     * To mark the read's mate as mapped, use {@link #setMatePosition}
     *
     * Calling this method has the additional effect of marking the read as paired, as if {@link #setIsPaired}
     * were invoked with true.
     */
    void setMateIsUnmapped();

    /**
     * @return True if this read is on the reverse strand as opposed to the forward strand, otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isReverseStrand();

    /**
     * Mark the read as being on the reverse (or forward) strand.
     *
     * @param isReverseStrand True if this read is on the reverse strand as opposed to the forward strand, otherwise false.
     */
    void setIsReverseStrand( final boolean isReverseStrand );

    /**
     * @return True if this read's mate is on the reverse strand as opposed to the forward strand, otherwise false.
     * @throws IllegalStateException if the read is not paired (has no mate)
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean mateIsReverseStrand();

    /**
     * Mark the read's mate as being on the reverse (or forward) strand.
     *
     * Calling this method has the additional effect of marking the read as paired, as if {@link #setIsPaired}
     * were invoked with true.
     *
     * @param mateIsReverseStrand True if this read's mate is on the reverse strand as opposed to the forward strand,
     *                            otherwise false.
     */
    void setMateIsReverseStrand( final boolean mateIsReverseStrand );

    /**
     * @return True if this read is paired and is the first read in the pair, otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isFirstOfPair();

    /**
     * Mark the read as the first read of a pair.
     *
     * Calling this method has the additional effects of marking the read as paired, as if {@link #setIsPaired}
     * were invoked with true, and also marks the read as NOT being the second of a pair.
     */
    void setIsFirstOfPair();

    /**
     * @return True if this read is paired and is the second read in the pair, otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isSecondOfPair();

    /**
     * Mark the read as the second read of a pair.
     *
     * Calling this method has the additional effects of marking the read as paired, as if {@link #setIsPaired}
     * were invoked with true, and also marks the read as NOT being the first of a pair.
     */
    void setIsSecondOfPair();

    /**
     * @return True if this is a secondary alignment (an alternative to the primary alignment), otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isSecondaryAlignment();

    /**
     * Mark the read as a secondary alignment (an alternative to the primary alignment)
     *
     * @param isSecondaryAlignment True if this is a secondary alignment, otherwise false.
     */
    void setIsSecondaryAlignment( final boolean isSecondaryAlignment );

    /**
     * @return True if this is a supplementary alignment (used in the representation of a chimeric alignment), otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isSupplementaryAlignment();

    /**
     * Mark the read as a supplementary alignment (used in the representation of a chimeric alignment)
     *
     * @param isSupplementaryAlignment True if this is a supplementary alignment, otherwise false.
     */
    void setIsSupplementaryAlignment( final boolean isSupplementaryAlignment );

    /**
     * @return True if this read fails platform/vendor quality checks, otherwise false
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean failsVendorQualityCheck();

    /**
     * Mark the read as failing platform/vendor quality checks
     *
     * @param failsVendorQualityCheck True if this read fails platform/vendor quality checks, otherwise false
     */
    void setFailsVendorQualityCheck( final boolean failsVendorQualityCheck );

    /**
     * @return True if this read is a PCR or optical duplicate, otherwise false.
     * @throws GATKException.MissingReadField if this information is not available
     */
    boolean isDuplicate();

    /**
     * Mark the read as a PCR or optical duplicate.
     *
     * @param isDuplicate True if this read is a PCR or optical duplicate, otherwise false.
     */
    void setIsDuplicate( final boolean isDuplicate );

    /**
     * Check whether this read has a particular attribute
     *
     * @param attributeName name of the attribute to search for
     * @return true if the read has an attribute with the given name, otherwise false
     */
    boolean hasAttribute( final String attributeName );

    /**
     * Retrieve the value of a particular attribute typed as an integer.
     *
     * @param attributeName name of the attribute to retrieve
     * @return integer value of the requested attribute, or {@code null} if the attribute is not present
     * @throws GATKException.ReadAttributeTypeMismatch if the attribute
     *         value cannot be typed as an integer
     */
    Integer getAttributeAsInteger( final String attributeName );

    /**
     * Retrieve the value of a particular attribute typed as a String.
     *
     * @param attributeName name of the attribute to retrieve
     * @return String value of the requested attribute, or {@code null} if the attribute is not present
     * @throws GATKException.ReadAttributeTypeMismatch if the attribute
     *         value cannot be typed as a single String value.
     */
    String getAttributeAsString( final String attributeName );

    /**
     * Retrieve the value of a particular attribute typed as a byte array.
     *
     * Makes a defensive copy of an existing byte array within the read if necessary, so modifying
     * the return value will not modify the attribute value within the read.
     *
     * @param attributeName name of the attribute to retrieve
     * @return byte array value of the requested attribute, or {@code null} if the attribute is not present
     * @throws GATKException.ReadAttributeTypeMismatch if the attribute
     *         value cannot be typed as a byte array.
     */
    byte[] getAttributeAsByteArray( final String attributeName );

    /**
     * Set an integer-valued attribute on the read.
     *
     * @param attributeName Name of the attribute to set. Must be legal according to {@link ReadUtils#assertAttributeNameIsLegal}
     * @param attributeValue Integer value of the attribute (may be {@code null})
     * @throws IllegalArgumentException if the attribute name is illegal according to {@link ReadUtils#assertAttributeNameIsLegal}
     */
    void setAttribute( final String attributeName, final Integer attributeValue );

    /**
     * Set a String-valued attribute on the read.
     *
     * @param attributeName Name of the attribute to set. Must be legal according to {@link ReadUtils#assertAttributeNameIsLegal}
     * @param attributeValue String value of the attribute (may be {@code null})
     * @throws IllegalArgumentException if the attribute name is illegal according to {@link ReadUtils#assertAttributeNameIsLegal}
     */
    void setAttribute( final String attributeName, final String attributeValue );

    /**
     * Set a byte array attribute on the read.
     *
     * @param attributeName Name of the attribute to set. Must be legal according to {@link ReadUtils#assertAttributeNameIsLegal}
     * @param attributeValue byte array value of the attribute (may be {@code null} or empty)
     * @throws IllegalArgumentException if the attribute name is illegal according to {@link ReadUtils#assertAttributeNameIsLegal}
     */
    void setAttribute( final String attributeName, final byte[] attributeValue );

    /**
     * Clear an individual attribute on the read.
     *
     * @param attributeName Name of the attribute to clear. Must be legal according to {@link ReadUtils#assertAttributeNameIsLegal}
     * @throws IllegalArgumentException if the attribute name is illegal according to {@link ReadUtils#assertAttributeNameIsLegal}
     */
    void clearAttribute( final String attributeName );

    /**
     * Clear all attributes on the read.
     */
    void clearAttributes();

    /**
     * @return A copy of this read. The copy will not necessarily be a true deep copy (the fields inside the
     *         encapsulated read itself may be shallow copied), but should be safe to use freely in general
     *         given that all GATKRead methods that return mutable reference types make defensive copies
     *         (with the exception of the conversion methods {@link #convertToSAMRecord} and
     *         {@link #convertToGoogleGenomicsRead}, but these are safe to call on copies since the
     *         encapsulated reads do get shallow copied at a minimum by this method, so modifications
     *         to the fields within a copied read will not alter the original).
     *
     * TODO: guarantee a deep copy here (requires writing a deep copy method for SAMRecord; copy() for
     * TODO: Google Genomics reads already performs a deep copy)
     */
    GATKRead copy();

    /**
     * Convert this read into a SAMRecord.
     *
     * Warning: the return value is not guaranteed to be independent from this read (eg., if the read
     * is already in SAMRecord format, no copy will be made).
     *
     * @param header required header for the SAMRecord
     * @return This read as a SAMRecord
     */
    SAMRecord convertToSAMRecord( final SAMFileHeader header );

    /**
     * Convert this read into a Google Genomics model read.
     *
     * Warning: the return value is not guaranteed to be independent from this read (eg., if the read
     * is already in Google Genomics format, no copy will be made).
     *
     * @return This read as a Google Genomics model read.
     */
    Read convertToGoogleGenomicsRead();

    /**
     * @param other Object to compare against
     * @return True if this is the same class as other and the underlying reads are equal, ignoring UUIDs in the comparison
     */
    boolean equalsIgnoreUUID( final Object other );

}

