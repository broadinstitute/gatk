package org.broadinstitute.hellbender.utils.read;


import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Implementation of the {@link GATKRead} interface for the {@link SAMRecord} class.
 *
 * This adapter wraps a {@link SAMRecord} without making a copy, so construction is cheap,
 * but care must be exercised if the underlying read has been exposed somewhere before
 * wrapping.
 */
public class SAMRecordToGATKReadAdapter implements GATKRead, Serializable {
    private static final long serialVersionUID = 1L;

    private final SAMRecord samRecord;

    public SAMRecordToGATKReadAdapter( final SAMRecord samRecord ) {
        this.samRecord = samRecord;
    }

    /**
     * Produces a SAMRecordToGATKReadAdapter wrapping the provided SAMRecord,
     * and nulls out the header in the encapsulated read. This is useful for
     * Spark tools, in order to avoid serializing the SAMFileHeader for each
     * record.
     *
     * @param samRecord read to adapt (header will be stripped)
     * @return SAMRecordToGATKReadAdapter wrapping the headerless samRecord
     */
    public static SAMRecordToGATKReadAdapter headerlessReadAdapter( final SAMRecord samRecord ) {
        samRecord.setHeaderStrict(null);
        return new SAMRecordToGATKReadAdapter(samRecord);
    }

    @Override
    public String getName() {
        return samRecord.getReadName();
    }

    @Override
    public void setName( final String name ) {
        samRecord.setReadName(name);
    }

    @Override
    public String getContig() {
        if ( isUnmapped() ) {
            return null;
        }

        // Guaranteed not to be null or SAMRecord.NO_ALIGNMENT_REFERENCE_NAME due to the isUnmapped() check above
        return samRecord.getReferenceName();
    }

    @Override
    public int getStart() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        // Guaranteed not to be SAMRecord.NO_ALIGNMENT_START due to the isUnmapped() check above
        return samRecord.getAlignmentStart();
    }

    @Override
    public int getEnd() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        // Guaranteed not to be SAMRecord.NO_ALIGNMENT_START due to the isUnmapped() check above
        return samRecord.getAlignmentEnd();
    }

    @Override
    public void setPosition( final String contig, final int start ) {
        if ( contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1 ) {
            throw new IllegalArgumentException("contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME + ", and start must be >= 1");
        }

        samRecord.setReferenceName(contig);
        samRecord.setAlignmentStart(start);
        samRecord.setReadUnmappedFlag(false);
    }

    @Override
    public void setPosition( final Locatable locatable ) {
        if ( locatable == null ) {
            throw new IllegalArgumentException("Cannot set read position to null");
        }

        setPosition(locatable.getContig(), locatable.getStart());
    }

    @Override
    public int getUnclippedStart() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        return samRecord.getUnclippedStart();
    }

    @Override
    public int getUnclippedEnd() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        return samRecord.getUnclippedEnd();
    }

    @Override
    public String getMateContig() {
        if ( mateIsUnmapped() ) {
            return null;
        }

        return samRecord.getMateReferenceName();
    }

    @Override
    public int getMateStart() {
        if ( mateIsUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        return samRecord.getMateAlignmentStart();
    }

    @Override
    public void setMatePosition( final String contig, final int start ) {
        if ( contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1 ) {
            throw new IllegalArgumentException("contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME + ", and start must be >= 1");
        }

        // Calling this method has the additional effect of marking the read as paired
        setIsPaired(true);

        samRecord.setMateReferenceName(contig);
        samRecord.setMateAlignmentStart(start);
        samRecord.setMateUnmappedFlag(false);
    }

    @Override
    public void setMatePosition( final Locatable locatable ) {
        Utils.nonNull(locatable, "Cannot set mate position to null");

        setMatePosition(locatable.getContig(), locatable.getStart());
    }

    @Override
    public int getFragmentLength() {
        return samRecord.getInferredInsertSize();
    }

    @Override
    public void setFragmentLength( final int fragmentLength ) {
        // May be negative if mate maps to lower position than read
        samRecord.setInferredInsertSize(fragmentLength);
    }

    @Override
    public int getMappingQuality() {
        return samRecord.getMappingQuality() != SAMRecord.NO_MAPPING_QUALITY ? samRecord.getMappingQuality() : ReadConstants.NO_MAPPING_QUALITY;
    }

    @Override
    public void setMappingQuality( final int mappingQuality ) {
        if ( mappingQuality < 0 || mappingQuality > 255 ) {
            throw new IllegalArgumentException("mapping quality must be >= 0 and <= 255");
        }

        samRecord.setMappingQuality(mappingQuality);
    }

    @Override
    public byte[] getBases() {
        final byte[] bases = samRecord.getReadBases();

        // Make a defensive copy to protect against direct modification of the returned array
        return bases != null ? Arrays.copyOf(bases, bases.length) : new byte[0];
    }

    //Overridden default method to avoid a call to getBases which makes a copy of data
    @Override
    public byte getBase(final int i){
        final byte[] bases = samRecord.getReadBases();
        if (bases == null){
            throw new IllegalArgumentException("Invalid call - there are no bases");
        }
        Utils.validIndex(i, bases.length);
        return bases[i];
    }

    @Override
    public int getLength() {
        final byte[] bases = samRecord.getReadBases();
        return bases == null ? 0 : bases.length;
    }

    @Override
    public void setBases( final byte[] bases ) {
        samRecord.setReadBases(bases);
    }

    @Override
    public byte[] getBaseQualities() {
        final byte[] baseQualities = samRecord.getBaseQualities();

        // Make a defensive copy to protect against direct modification of the returned array
        return baseQualities != null ? Arrays.copyOf(baseQualities, baseQualities.length) : new byte[0];
    }

    @Override
    public int getBaseQualityCount(){
        final byte[] baseQualities = samRecord.getBaseQualities();
        return baseQualities == null ? 0 : baseQualities.length;
    }

    //Overridden default method to avoid a call to getBaseQualities which makes a copy of data
    @Override
    public byte getBaseQuality(final int i){
        final byte[] baseQualities = samRecord.getBaseQualities();
        if (baseQualities == null){
            throw new IllegalArgumentException("Invalid call - there are no baseQualities");
        }
        Utils.validIndex(i, baseQualities.length);
        return baseQualities[i];
    }

    @Override
    public void setBaseQualities( final byte[] baseQualities ) {
        if ( baseQualities != null ) {
            for ( byte b : baseQualities ) {
                if ( b < 0 ) {
                    throw new IllegalArgumentException("Base quality score " + b + " is invalid");
                }
            }
        }

        samRecord.setBaseQualities(baseQualities);
    }

    @Override
    public Cigar getCigar() {
        // Make a defensive copy before returning to guard against modification of the return value,
        // since Cigar is a mutable type:
        return samRecord.getCigar() != null ? new Cigar(samRecord.getCigar().getCigarElements()) : new Cigar();
    }

    /**
     * This implementation does not make a new Cigar object but instead provides
     * an unmodifiable view of the underlying list of CigarElements.
     * This is done to reduce the amount of object allocation.
     */
    @Override
    public List<CigarElement> getCigarElements(){
        //Cigar.getCigarElements returns an unmodifiable list so we don't wrap it again
        return samRecord.getCigar() == null ? Collections.emptyList() : samRecord.getCigar().getCigarElements();
    }

    /**
     * This implementation saves time by not making an unmodifiable view of the list of
     * elements but returns the length of the list directly (or 0 if there's no cigar).
     */
    @Override
    public int numCigarElements(){
        return samRecord.getCigar() == null ? 0 : samRecord.getCigar().numCigarElements();
    }

    @Override
    public void setCigar( final Cigar cigar ) {
        samRecord.setCigar(cigar);
    }

    @Override
    public void setCigar( final String cigarString ) {
        samRecord.setCigarString(cigarString);
    }

    @Override
    public String getReadGroup() {
        // May return null
        return (String)samRecord.getAttribute(SAMTagUtil.getSingleton().RG);
    }

    @Override
    public void setReadGroup( final String readGroupID ) {
        samRecord.setAttribute(SAMTag.RG.name(), readGroupID);
    }

    @Override
    public boolean isPaired() {
        return samRecord.getReadPairedFlag();
    }

    @Override
    public void setIsPaired( final boolean isPaired ) {
        samRecord.setReadPairedFlag(isPaired);
        if ( ! isPaired ) {
            samRecord.setProperPairFlag(false);
        }
    }

    @Override
    public boolean isProperlyPaired() {
        return isPaired() && samRecord.getProperPairFlag();
    }

    @Override
    public void setIsProperlyPaired( final boolean isProperlyPaired ) {
        if ( isProperlyPaired ) {
            setIsPaired(true);
        }

        samRecord.setProperPairFlag(isProperlyPaired);
    }

    @Override
    public boolean isUnmapped() {
        return samRecord.getReadUnmappedFlag() ||
               samRecord.getReferenceName() == null || samRecord.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
               samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    @Override
    public void setIsUnmapped() {
        samRecord.setReadUnmappedFlag(true);
    }

    @Override
    public boolean mateIsUnmapped() {
        if ( ! isPaired() ) {
            throw new IllegalStateException("Cannot get mate information for an unpaired read");
        }

        return samRecord.getMateUnmappedFlag() ||
               samRecord.getMateReferenceName() == null || samRecord.getMateReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
               samRecord.getMateAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    @Override
    public void setMateIsUnmapped() {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setMateUnmappedFlag(true);
    }

    @Override
    public boolean isReverseStrand() {
        return samRecord.getReadNegativeStrandFlag();
    }

    @Override
    public void setIsReverseStrand( final boolean isReverseStrand ) {
        samRecord.setReadNegativeStrandFlag(isReverseStrand);
    }

    @Override
    public boolean mateIsReverseStrand() {
        if ( ! isPaired() ) {
            throw new IllegalStateException("Cannot get mate information for an unpaired read");
        }

        return samRecord.getMateNegativeStrandFlag();
    }

    @Override
    public void setMateIsReverseStrand( final boolean mateIsReverseStrand ) {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setMateNegativeStrandFlag(mateIsReverseStrand);
    }

    @Override
    public boolean isFirstOfPair() {
        return isPaired() && samRecord.getFirstOfPairFlag();
    }

    @Override
    public void setIsFirstOfPair() {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setFirstOfPairFlag(true);
        samRecord.setSecondOfPairFlag(false);
    }

    @Override
    public boolean isSecondOfPair() {
        return isPaired() && samRecord.getSecondOfPairFlag();
    }

    @Override
    public void setIsSecondOfPair() {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        samRecord.setSecondOfPairFlag(true);
        samRecord.setFirstOfPairFlag(false);
    }

    @Override
    public boolean isSecondaryAlignment() {
        return samRecord.getNotPrimaryAlignmentFlag();
    }

    @Override
    public void setIsSecondaryAlignment( final boolean isSecondaryAlignment ) {
        samRecord.setNotPrimaryAlignmentFlag(isSecondaryAlignment);
    }

    @Override
    public boolean isSupplementaryAlignment() {
        return samRecord.getSupplementaryAlignmentFlag();
    }

    @Override
    public void setIsSupplementaryAlignment( final boolean isSupplementaryAlignment ) {
        samRecord.setSupplementaryAlignmentFlag(isSupplementaryAlignment);
    }

    @Override
    public boolean failsVendorQualityCheck() {
        return samRecord.getReadFailsVendorQualityCheckFlag();
    }

    @Override
    public void setFailsVendorQualityCheck( final boolean failsVendorQualityCheck ) {
        samRecord.setReadFailsVendorQualityCheckFlag(failsVendorQualityCheck);
    }

    @Override
    public boolean isDuplicate() {
        return samRecord.getDuplicateReadFlag();
    }

    @Override
    public void setIsDuplicate( final boolean isDuplicate ) {
        samRecord.setDuplicateReadFlag(isDuplicate);
    }

    @Override
    public boolean hasAttribute( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        return samRecord.getAttribute(attributeName) != null;
    }

    @Override
    public Integer getAttributeAsInteger( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = samRecord.getAttribute(attributeName);

        if ( attributeValue == null ) {
            return null;
        }
        else if ( attributeValue instanceof Integer ) {
            return (Integer)attributeValue;
        }
        else {
            try {
                return Integer.parseInt(attributeValue.toString());
            }
            catch ( NumberFormatException e ) {
                throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer", e);
            }
        }
    }

    @Override
    public String getAttributeAsString( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = samRecord.getAttribute(attributeName);

        return attributeValue != null ? attributeValue.toString() : null;
    }

    @Override
    public byte[] getAttributeAsByteArray( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = samRecord.getAttribute(attributeName);

        if ( attributeValue == null ) {
            return null;
        }
        else if ( attributeValue instanceof byte[] ) {
            // In the case where the attribute value is already a byte[], make a defensive
            // copy before returning to guard against modification of the return value.
            final byte[] ret = (byte[])attributeValue;
            return Arrays.copyOf(ret, ret.length);
        }
        else if ( attributeValue instanceof String ) {
            return ((String)attributeValue).getBytes();
        }
        else {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, "byte array");
        }
    }

    @Override
    public void setAttribute( final String attributeName, final Integer attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        samRecord.setAttribute(attributeName, attributeValue);
    }

    @Override
    public void setAttribute( final String attributeName, final String attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        samRecord.setAttribute(attributeName, attributeValue);
    }

    @Override
    public void setAttribute( final String attributeName, final byte[] attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        samRecord.setAttribute(attributeName, attributeValue);
    }

    @Override
    public void clearAttribute( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        samRecord.setAttribute(attributeName, null);
    }

    @Override
    public void clearAttributes() {
        samRecord.clearAttributes();
    }

    @Override
    public GATKRead copy() {
        // Produces a shallow but "safe to use" copy.
        return new SAMRecordToGATKReadAdapter(ReadUtils.cloneSAMRecord(samRecord));
    }

    @Override
    public GATKRead deepCopy() {
        // Produces a true deep copy.
        return new SAMRecordToGATKReadAdapter(samRecord.deepCopy());
    }

    @Override
    public String getSAMString() {
        return samRecord.getSAMString();
    }

    @Override
    public SAMRecord convertToSAMRecord( final SAMFileHeader header ) {
        samRecord.setHeaderStrict(header);
        return samRecord;
    }

    public SAMRecord getEncapsulatedSamRecord() {
        return samRecord;
    }

    @Override
    public Read convertToGoogleGenomicsRead() {
        // TODO: this converter is imperfect/lossy and should either be patched or replaced
        return com.google.cloud.genomics.utils.ReadUtils.makeRead(samRecord);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SAMRecordToGATKReadAdapter that = (SAMRecordToGATKReadAdapter) o;

        return !(samRecord != null ? !samRecord.equals(that.samRecord) : that.samRecord != null);

    }

    @Override
    public int hashCode() {
        return samRecord != null ? samRecord.hashCode() : 0;
    }

    @Override
    public String toString() {
        return commonToString();
    }

    public boolean hasHeader() {
        return samRecord.getHeader() != null;
    }

    public void setHeader(SAMFileHeader header) {
        samRecord.setHeaderStrict(header);
    }
}
