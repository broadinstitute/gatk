package org.broadinstitute.hellbender.utils.read;


import com.google.api.services.genomics.model.LinearAlignment;
import com.google.api.services.genomics.model.Position;
import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.Serializable;
import java.util.*;

/**
 * Implementation of the {@link GATKRead} interface for the Google Genomics {@link Read} class.
 *
 * This adapter wraps a {@link Read} without making a copy, so construction is cheap,
 * but care must be exercised if the underlying read has been exposed somewhere before
 * wrapping.
 */
public final class GoogleGenomicsReadToGATKReadAdapter implements GATKRead, Serializable {
    private static final long serialVersionUID = 1L;

    private final Read genomicsRead;

    private final static String SAM_FIELD_SEPARATOR = "\t";
    private final static String SAM_ATTRIBUTE_SEPARATOR = ":";

    public GoogleGenomicsReadToGATKReadAdapter( final Read genomicsRead ) {
        this.genomicsRead = genomicsRead;
    }

    private static <T> T assertFieldValueNotNull( final T fieldValue, final String fieldName ) {
        if ( fieldValue == null ) {
            throw new GATKException.MissingReadField(fieldName);
        }
        return fieldValue;
    }

    private void assertHasAlignment() {
        assertFieldValueNotNull(genomicsRead.getAlignment(), "alignment");
    }

    private void assertHasPosition() {
        assertHasAlignment();
        assertFieldValueNotNull(genomicsRead.getAlignment().getPosition(), "position");
    }

    private boolean positionIsUnmapped( final Position position ) {
        return position == null ||
               position.getReferenceName() == null || position.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
               position.getPosition() == null || position.getPosition() < 0;
    }

    private void makeAlignmentIfNecessary() {
        if ( genomicsRead.getAlignment() == null ) {
            genomicsRead.setAlignment(new LinearAlignment());
        }
    }

    private void makePositionIfNecessary() {
        makeAlignmentIfNecessary();

        if ( genomicsRead.getAlignment().getPosition() == null ) {
            genomicsRead.getAlignment().setPosition(new Position());
        }
    }

    private void makeMatePositionIfNecessary() {
        if ( genomicsRead.getNextMatePosition() == null ) {
            genomicsRead.setNextMatePosition(new Position());
        }
    }

    @Override
    public String getName() {
        return genomicsRead.getFragmentName();
    }

    @Override
    public void setName( final String name ) {
        genomicsRead.setFragmentName(name);
    }

    @Override
    public String getContig() {
        if ( isUnmapped() ) {
            return null;
        }

        // Guaranteed non-null due to isUnmapped() check above.
        return genomicsRead.getAlignment().getPosition().getReferenceName();
    }

    @Override
    public int getStart() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        // Guaranteed non-null due to isUnmapped() check above.
        // Convert from 0-based to 1-based start position
        return genomicsRead.getAlignment().getPosition().getPosition().intValue() + 1;
    }

    @Override
    public int getEnd() {
        if ( isUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        // Guaranteed non-null due to isUnmapped() check above.
        // Position in genomicsRead is 0-based, so add getCigar().getReferenceLength() to it,
        // not getCigar().getReferenceLength() - 1, in order to convert to a 1-based end position.
        return genomicsRead.getAlignment().getPosition().getPosition().intValue() + getCigar().getReferenceLength();
    }

    @Override
    public void setPosition( final String contig, final int start ) {
        if ( contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1 ) {
            throw new IllegalArgumentException("contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME + ", and start must be >= 1");
        }

        // NOTE: this implementation will retain the previous value of the mate reverseStrandFlag if one
        // is already set. We should probably unset it here. However, once it is unset, an exception will
        // be thrown when querying the value, and I don't think we have a guard method that can tell if its
        // safe to call.
        makePositionIfNecessary();

        genomicsRead.getAlignment().getPosition().setReferenceName(contig);
        // Convert from a 1-based to a 0-based position
        genomicsRead.getAlignment().getPosition().setPosition((long)start - 1);
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
        final int start = getStart();
        return start == ReadConstants.UNSET_POSITION ? ReadConstants.UNSET_POSITION :
                                                       SAMUtils.getUnclippedStart(start, getCigar());
    }

    @Override
    public int getUnclippedEnd() {
        final int end = getEnd();
        return end == ReadConstants.UNSET_POSITION ? ReadConstants.UNSET_POSITION :
                SAMUtils.getUnclippedEnd(getEnd(), getCigar());
    }

    @Override
    public String getMateContig() {
        if ( mateIsUnmapped() ) {
            return null;
        }

        // Guaranteed non-null due to mateIsUnmapped() check above.
        return genomicsRead.getNextMatePosition().getReferenceName();
    }

    @Override
    public int getMateStart() {
        if ( mateIsUnmapped() ) {
            return ReadConstants.UNSET_POSITION;
        }

        // Guaranteed non-null due to mateIsUnmapped() check above.
        // Convert from 0-based to 1-based position.
        return genomicsRead.getNextMatePosition().getPosition().intValue() + 1;
    }

    @Override
    public void setMatePosition( final String contig, final int start ) {
        if ( contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1 ) {
            throw new IllegalArgumentException("contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME + ", and start must be >= 1");
        }

        // Calling this method has the additional effect of marking the read as paired
        setIsPaired(true);

        // NOTE: this implementation will retain the previous value of the mate reverseStrandFlag if one
        // is already set. We should probably unset it here. However, once it is unset, an exception will
        // be thrown when querying the value, and I don't think we have a guard method that can tell if its
        // safe to call.
        makeMatePositionIfNecessary();
        genomicsRead.getNextMatePosition().setReferenceName(contig);
        // Convert from a 1-based to a 0-based position
        genomicsRead.getNextMatePosition().setPosition((long)start - 1);
    }

    @Override
    public void setMatePosition( final Locatable locatable ) {
        if ( locatable == null ) {
            throw new IllegalArgumentException("Cannot set mate position to null");
        }

        setMatePosition(locatable.getContig(), locatable.getStart());
    }

    @Override
    public int getFragmentLength() {
        return genomicsRead.getFragmentLength() != null ? genomicsRead.getFragmentLength() : 0;
    }

    @Override
    public void setFragmentLength( final int fragmentLength ) {
        // May be negative if mate maps to lower position than read
        genomicsRead.setFragmentLength(fragmentLength);
    }

    @Override
    public int getMappingQuality() {
        if ( genomicsRead.getAlignment() == null || genomicsRead.getAlignment().getMappingQuality() == null ) {
            return ReadConstants.NO_MAPPING_QUALITY;
        }

        return genomicsRead.getAlignment().getMappingQuality();
    }

    @Override
    public void setMappingQuality( final int mappingQuality ) {
        if ( mappingQuality < 0 || mappingQuality > 255 ) {
            throw new IllegalArgumentException("mapping quality must be >= 0 and <= 255");
        }

        makeAlignmentIfNecessary();
        genomicsRead.getAlignment().setMappingQuality(mappingQuality);
    }

    @Override
    public byte[] getBases() {
        final String basesString = genomicsRead.getAlignedSequence();
        if ( getLength() == 0) {
            return new byte[0];
        }

        return StringUtil.stringToBytes(basesString);
    }

    @Override
    public int getLength() {
        final String basesString = genomicsRead.getAlignedSequence();
        if (basesString == null || basesString.equals(SAMRecord.NULL_SEQUENCE_STRING)){
            return 0;
        } else {
            return basesString.length();
        }
    }

    @Override
    public void setBases( final byte[] bases ) {
        genomicsRead.setAlignedSequence(bases != null ? StringUtil.bytesToString(bases) : null);
    }

    @Override
    public byte[] getBaseQualities() {
        final List<Integer> baseQualities = genomicsRead.getAlignedQuality();
        if ( baseQualities == null || baseQualities.isEmpty() ) {
            return new byte[0];
        }

        byte[] convertedBaseQualities = new byte[baseQualities.size()];
        for ( int i = 0; i < baseQualities.size(); ++i ) {
            if ( baseQualities.get(i) < 0 || baseQualities.get(i) > Byte.MAX_VALUE ) {
                throw new GATKException("Base quality score " + baseQualities.get(i) + " is invalid and/or not convertible to byte");
            }
            convertedBaseQualities[i] = baseQualities.get(i).byteValue();
        }

        return convertedBaseQualities;
    }

    @Override
    public void setBaseQualities( final byte[] baseQualities ) {
        if ( baseQualities == null ) {
            genomicsRead.setAlignedQuality(null);
            return;
        }

        final List<Integer> convertedBaseQualities = new ArrayList<>(baseQualities.length);
        for ( byte b : baseQualities ) {
            if ( b < 0 ) {
                throw new IllegalArgumentException("Base quality score " + b + " is invalid");
            }

            convertedBaseQualities.add((int)b);
        }

        genomicsRead.setAlignedQuality(convertedBaseQualities.isEmpty() ? null : convertedBaseQualities);
    }

    @Override
    public Cigar getCigar() {
        if ( genomicsRead.getAlignment() == null || genomicsRead.getAlignment().getCigar() == null ) {
            return new Cigar();
        }

        return CigarConversionUtils.convertCigarUnitListToSAMCigar(genomicsRead.getAlignment().getCigar());
    }

    @Override
    public void setCigar( final Cigar cigar ) {
        makeAlignmentIfNecessary();
        genomicsRead.getAlignment().setCigar(cigar != null ? CigarConversionUtils.convertSAMCigarToCigarUnitList(cigar) : null);
    }

    @Override
    public void setCigar( final String cigarString ) {
        makeAlignmentIfNecessary();
        genomicsRead.getAlignment().setCigar(cigarString != null ? CigarConversionUtils.convertSAMCigarToCigarUnitList(TextCigarCodec.decode(cigarString)) : null);
    }

    @Override
    public String getReadGroup() {
        return genomicsRead.getReadGroupId();
    }

    @Override
    public void setReadGroup( final String readGroupID ) {
        genomicsRead.setReadGroupId(readGroupID);
    }

    @Override
    public boolean isPaired() {
        assertFieldValueNotNull(genomicsRead.getNumberReads(), "number of reads");
        return genomicsRead.getNumberReads() == 2;
    }

    @Override
    public void setIsPaired( final boolean isPaired ) {
        // Setting isPaired==false changes numberOfReads to 1, but the GenomicsConverter will
        // still set the isFirstOfPair flag unless the readNumber is 1.
        // https://github.com/broadinstitute/hellbender/issues/891
        genomicsRead.setNumberReads(isPaired ? 2 : 1);
        if ( ! isPaired ) {
            genomicsRead.setProperPlacement(false);
        }
    }

    @Override
    public boolean isProperlyPaired() {
        assertFieldValueNotNull(genomicsRead.getProperPlacement(), "proper placement");
        return isPaired() && genomicsRead.getProperPlacement();
    }

    @Override
    public void setIsProperlyPaired( final boolean isProperlyPaired ) {
        if ( isProperlyPaired ) {
            setIsPaired(true);
        }

        genomicsRead.setProperPlacement(isProperlyPaired);
    }

    @Override
    public boolean isUnmapped() {
        return genomicsRead.getAlignment() == null ||
               positionIsUnmapped(genomicsRead.getAlignment().getPosition());
    }

    @Override
    public void setIsUnmapped() {
        genomicsRead.setAlignment(null);
    }

    @Override
    public boolean mateIsUnmapped() {
        if ( ! isPaired() ) {
            throw new IllegalStateException("Cannot get mate information for an unpaired read");
        }

        return positionIsUnmapped(genomicsRead.getNextMatePosition());
    }

    @Override
    public void setMateIsUnmapped() {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        genomicsRead.setNextMatePosition(null);
    }

    @Override
    public boolean isReverseStrand() {
        assertHasPosition();
        return assertFieldValueNotNull(genomicsRead.getAlignment().getPosition().getReverseStrand(), "strand");
    }

    @Override
    public void setIsReverseStrand( final boolean isReverseStrand ) {
        makePositionIfNecessary();
        genomicsRead.getAlignment().getPosition().setReverseStrand(isReverseStrand);
    }

    @Override
    public boolean mateIsReverseStrand() {
        if ( ! isPaired() ) {
            throw new IllegalStateException("Cannot get mate information for an unpaired read");
        }

        final Position matePosition = assertFieldValueNotNull(genomicsRead.getNextMatePosition(), "mate position");
        return assertFieldValueNotNull(matePosition.getReverseStrand(), "mate strand");
    }

    @Override
    public void setMateIsReverseStrand( final boolean mateIsReverseStrand ) {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        makeMatePositionIfNecessary();
        genomicsRead.getNextMatePosition().setReverseStrand(mateIsReverseStrand);
    }

    @Override
    public boolean isFirstOfPair() {
        final int readNumber = assertFieldValueNotNull(genomicsRead.getReadNumber(), "read number");
        return isPaired() && readNumber == 0;
    }

    @Override
    public void setIsFirstOfPair() {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        genomicsRead.setReadNumber(0);
    }

    @Override
    public boolean isSecondOfPair() {
        final int readNumber = assertFieldValueNotNull(genomicsRead.getReadNumber(), "read number");
        return isPaired() && readNumber == 1;
    }

    @Override
    public void setIsSecondOfPair() {
        // Calling this method has the side effect of marking the read as paired.
        setIsPaired(true);

        genomicsRead.setReadNumber(1);
    }

    @Override
    public boolean isSecondaryAlignment() {
        return assertFieldValueNotNull(genomicsRead.getSecondaryAlignment(), "secondary alignment");
    }

    @Override
    public void setIsSecondaryAlignment( final boolean isSecondaryAlignment ) {
        genomicsRead.setSecondaryAlignment(isSecondaryAlignment);
    }

    @Override
    public boolean isSupplementaryAlignment() {
        return assertFieldValueNotNull(genomicsRead.getSupplementaryAlignment(), "supplementary alignment");
    }

    @Override
    public void setIsSupplementaryAlignment( final boolean isSupplementaryAlignment ) {
        genomicsRead.setSupplementaryAlignment(isSupplementaryAlignment);
    }

    @Override
    public boolean failsVendorQualityCheck() {
        return assertFieldValueNotNull(genomicsRead.getFailedVendorQualityChecks(), "failed vendor quality checks");
    }

    @Override
    public void setFailsVendorQualityCheck( final boolean failsVendorQualityCheck ) {
        genomicsRead.setFailedVendorQualityChecks(failsVendorQualityCheck);
    }

    @Override
    public boolean isDuplicate() {
        return assertFieldValueNotNull(genomicsRead.getDuplicateFragment(), "duplicate fragment");
    }

    @Override
    public void setIsDuplicate( final boolean isDuplicate ) {
        genomicsRead.setDuplicateFragment(isDuplicate);
    }

    @Override
    public boolean hasAttribute( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        return genomicsRead.getInfo() != null && genomicsRead.getInfo().containsKey(attributeName);
    }

    private String getRawAttributeValue( final String attributeName, final String targetType ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);

        if ( genomicsRead.getInfo() == null ) {
            return null;
        }

        final List<String> rawValue = genomicsRead.getInfo().get(attributeName);
        if ( rawValue == null || rawValue.isEmpty() || rawValue.get(0) == null ) {
            return null;
        }

        // We don't support decoding attribute values represented as multiple Strings
        if ( rawValue.size() > 1 ) {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, targetType);
        }

        return rawValue.get(0);
    }

    @Override
    public Integer getAttributeAsInteger( final String attributeName ) {
        try {
            // Assume that integer attributes are encoded as a single String in the first position of the List of values for an attribute
            final String rawValue = getRawAttributeValue(attributeName, "integer");
            return rawValue != null ? Integer.parseInt(rawValue) : null;
        }
        catch ( NumberFormatException e ) {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer", e);
        }
    }

    @Override
    public String getAttributeAsString( final String attributeName ) {
        // Assume that String attributes are encoded as a single String in the first position of the List of values for an attribute
        return getRawAttributeValue(attributeName, "String");
    }

    @Override
    public byte[] getAttributeAsByteArray( final String attributeName ) {
        // Assume that byte array attributes are encoded as a single String in the first position of the List of values for an attribute,
        // and that the bytes of this String are directly convertible to byte array.
        final String rawValue = getRawAttributeValue(attributeName, "byte array");
        return rawValue != null ? rawValue.getBytes() : null;
    }

    private void makeInfoMapIfNecessary() {
        if ( genomicsRead.getInfo() == null ) {
            genomicsRead.setInfo(new LinkedHashMap<>());
        }
    }

    @Override
    public void setAttribute( final String attributeName, final Integer attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        makeInfoMapIfNecessary();
        if ( attributeValue == null ) {
            clearAttribute(attributeName);
            return;
        }

        final List<String> encodedValue = Arrays.asList(attributeValue.toString());
        genomicsRead.getInfo().put(attributeName, encodedValue);
    }

    @Override
    public void setAttribute( final String attributeName, final String attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        makeInfoMapIfNecessary();
        if ( attributeValue == null ) {
            clearAttribute(attributeName);
            return;
        }

        final List<String> encodedValue = Arrays.asList(attributeValue);
        genomicsRead.getInfo().put(attributeName, encodedValue);
    }

    @Override
    public void setAttribute( final String attributeName, final byte[] attributeValue ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        makeInfoMapIfNecessary();
        if ( attributeValue == null ) {
            clearAttribute(attributeName);
            return;
        }

        final List<String> encodedValue = Arrays.asList(new String(attributeValue));
        genomicsRead.getInfo().put(attributeName, encodedValue);
    }

    @Override
    public void clearAttribute( final String attributeName ) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);

        if ( genomicsRead.getInfo() != null ) {
            genomicsRead.getInfo().remove(attributeName);
        }
    }

    @Override
    public void clearAttributes() {
        genomicsRead.setInfo(null);
    }

    @Override
    public GATKRead copy() {

        workAroundHttpClientCloneBug();

        // clone() actually makes a deep copy of all fields here (via GenericData.clone())
        return new GoogleGenomicsReadToGATKReadAdapter(genomicsRead.clone());
    }

    @Override
    public GATKRead deepCopy() {
        // Produces a deep copy..
        return copy();
    }

    /**
     * Return a SAM string representation of this read.
     */
    @Override
    public String getSAMString() {
        final String contigName = getContig();

        return new StringBuilder().append(getName())
                .append(SAM_FIELD_SEPARATOR)
                .append(Integer.toString(ReadUtils.getSAMFlagsForRead(this)))
                .append(SAM_FIELD_SEPARATOR)
                .append(null == contigName ? SAMRecord.NO_ALIGNMENT_REFERENCE_NAME : contigName)
                .append(SAM_FIELD_SEPARATOR)
                .append(Integer.toString(getStart()))
                .append(SAM_FIELD_SEPARATOR)
                .append(Integer.toString(getMappingQuality()))
                .append(SAM_FIELD_SEPARATOR)
                .append(getCigar().toString())
                .append(SAM_FIELD_SEPARATOR)
                .append(getMateContigDisplayString(contigName))
                .append(SAM_FIELD_SEPARATOR)
                .append(isPaired() ? Integer.toString(getMateStart()) : ReadConstants.UNSET_POSITION)
                .append(SAM_FIELD_SEPARATOR)
                .append(Integer.toString(getFragmentLength())) // getInferredInsertSize
                .append(SAM_FIELD_SEPARATOR)
                .append(getBasesString())
                .append(SAM_FIELD_SEPARATOR)
                .append(ReadUtils.getBaseQualityString(this))
                .append(getAttributesDisplayString()).append("\n").toString();
    }

    @Override
    public SAMRecord convertToSAMRecord( final SAMFileHeader header ) {
        // TODO: this converter is imperfect and should either be patched or replaced completely.
        return GenomicsConverter.makeSAMRecord(genomicsRead, header);
    }

    @Override
    public Read convertToGoogleGenomicsRead() {
        return genomicsRead;
    }

    @Override
    public boolean equals( Object other ) {
        if ( this == other ) return true;
        if ( other == null || getClass() != other.getClass() ) return false;

        GoogleGenomicsReadToGATKReadAdapter that = (GoogleGenomicsReadToGATKReadAdapter)other;

        // The Read class does have a working equals() method (inherited from AbstractMap)
        return genomicsRead != null ? genomicsRead.equals(that.genomicsRead) : that.genomicsRead == null;
    }

    @Override
    public int hashCode() {
        return genomicsRead != null ? genomicsRead.hashCode() : 0;
    }

    private String getMateContigDisplayString(final String contigName) {
        String mateContigString;
        if (!isPaired()) {
            mateContigString = SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
        }
        else {
            final String mateContigName = getMateContig();
            if (null == mateContigName) {
                mateContigString = SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
            }
            else if (null != contigName &&
                    !contigName.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) &&
                    contigName.equals(mateContigName)) {
                mateContigString = "=";
            }
            else {
                mateContigString = mateContigName;
            }
        }
        return mateContigString;
    }

    private String getAttributesDisplayString() {
        final StringBuilder sb = new StringBuilder();
        if (null != genomicsRead.getInfo()) {
            final TreeSet<String> sortedAttributes = new TreeSet<>();
            final Set<String> keys = genomicsRead.getInfo().keySet();
            sortedAttributes.addAll(keys);
            // Add the RG attribute if necessary since its not returned by read.getInfo
            if (getReadGroup() != null) {
                sortedAttributes.add(SAMTag.RG.name());
            }
            for (String key: sortedAttributes) {
                sb.append(SAM_FIELD_SEPARATOR);
                sb.append(key
                        + SAM_ATTRIBUTE_SEPARATOR
                        + GenomicsConverter.getTagType(key)
                        + SAM_ATTRIBUTE_SEPARATOR
                        + resolveAttributeValue(keys, key)); // resolve readGroup value if necessary
            }
        }
        return sb.toString();
    }

    // TODO: Remove once https://github.com/broadinstitute/hellbender/issues/650 is solved.
    private void workAroundHttpClientCloneBug() {
        final List<Integer> alignedQuality = genomicsRead.getAlignedQuality();
        if (null!=alignedQuality) {
            genomicsRead.setAlignedQuality(new ArrayList<>(alignedQuality));
        }
        if (null!=genomicsRead.getInfo()) {
            Map<String, List<String>> infoCopy = new HashMap<>();
            for (Map.Entry<String, List<String>> entry : genomicsRead.getInfo().entrySet()) {
                infoCopy.put(entry.getKey(), new ArrayList<>(entry.getValue()));
            }
            genomicsRead.setInfo(infoCopy);
        }
        final LinearAlignment alignment = genomicsRead.getAlignment();
        if (null!=alignment) {
            alignment.setCigar(new ArrayList<>(alignment.getCigar()));
        }
    }

    private String resolveAttributeValue(Set<String> keys, String key) {
        String attributeValue;
        if (keys.contains(key)) {
            attributeValue = getAttributeAsString(key);
        } else if (key.equals(SAMTag.RG.name())) {
            attributeValue = getReadGroup();
        }
        else {
            throw new GATKException.ShouldNeverReachHereException("Unresolvable attribute value: " + key);
        }
        return attributeValue;
    }

    @Override
    public String toString() {
        return commonToString();
    }
}
