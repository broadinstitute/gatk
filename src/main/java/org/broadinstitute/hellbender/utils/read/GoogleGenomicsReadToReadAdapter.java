package org.broadinstitute.hellbender.utils.read;


import com.google.api.services.genomics.model.LinearAlignment;
import com.google.api.services.genomics.model.Position;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.util.ArrayList;
import java.util.List;

public final class GoogleGenomicsReadToReadAdapter implements MutableRead {

    private final com.google.api.services.genomics.model.Read genomicsRead;

    public GoogleGenomicsReadToReadAdapter( final com.google.api.services.genomics.model.Read genomicsRead ) {
        this.genomicsRead = genomicsRead;
    }

    @Override
    public String getContig() {
        if ( isUnmapped() ) {
            return SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
        }

        // Guaranteed non-null due to isUnmapped() check above
        return genomicsRead.getAlignment().getPosition().getReferenceName();
    }

    @Override
    public int getStart() {
        if ( isUnmapped() ) {
            return GenomicIndexUtil.UNSET_GENOMIC_LOCATION;
        }

        // Guaranteed non-null due to isUnmapped() check above
        // Convert from 0-based to 1-based start position
        return genomicsRead.getAlignment().getPosition().getPosition().intValue() + 1;
    }

    @Override
    public int getEnd() {
        if ( isUnmapped() ) {
            return GenomicIndexUtil.UNSET_GENOMIC_LOCATION;
        }

        // Guaranteed non-null due to isUnmapped() check above
        // Position in genomicsRead is 0-based, so add getCigar().getReferenceLength() to it,
        // not getCigar().getReferenceLength() - 1, in order to convert to a 1-based end position.
        return genomicsRead.getAlignment().getPosition().getPosition().intValue() + getCigar().getReferenceLength();
    }

    @Override
    public String getName() {
        return genomicsRead.getFragmentName();
    }

    @Override
    public int getLength() {
        return getBases().length;
    }

    @Override
    public byte[] getBases() {
        final String basesString = genomicsRead.getAlignedSequence();
        if ( basesString == null || basesString.isEmpty() || basesString.equals(SAMRecord.NULL_SEQUENCE_STRING) ) {
            return new byte[0];
        }

        final byte[] bases = StringUtil.stringToBytes(basesString);
        ReadUtils.normalizeBases(bases);
        return bases;
    }

    @Override
    public byte[] getBaseQualities() {
        final List<Integer> baseQualities = genomicsRead.getAlignedQuality();
        if ( baseQualities == null ) {
            return new byte[0];
        }

        final byte[] convertedBaseQualities = new byte[baseQualities.size()];
        for ( int i = 0; i < baseQualities.size(); ++i ) {
            if ( baseQualities.get(i) < 0 || baseQualities.get(i) > Byte.MAX_VALUE ) {
                throw new GATKException("Base quality score " + baseQualities.get(i) + " is invalid and/or not convertible to byte");
            }
            convertedBaseQualities[i] = baseQualities.get(i).byteValue();
        }
        return convertedBaseQualities;
    }

    @Override
    public Cigar getCigar() {
        if ( genomicsRead.getAlignment() == null || genomicsRead.getAlignment().getCigar() == null ) {
            return new Cigar();
        }

        return CigarConverter.convertCigarUnitListToSAMCigar(genomicsRead.getAlignment().getCigar());
    }

    @Override
    public boolean isUnmapped() {
        return genomicsRead.getAlignment() == null ||
               positionIsUnmapped(genomicsRead.getAlignment().getPosition());
    }

    private boolean positionIsUnmapped( final Position position ) {
        return position == null ||
               position.getReferenceName() == null || position.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
               position.getPosition() == null || position.getPosition() < 0;
    }

    @Override
    public void setName( String name ) {
        genomicsRead.setFragmentName(name);
    }

    @Override
    public void setBases( byte[] bases ) {
        genomicsRead.setAlignedSequence(StringUtil.bytesToString(bases));
    }

    @Override
    public void setBaseQualities( byte[] baseQualities ) {
        final List<Integer> convertedBaseQualities = new ArrayList<>(baseQualities.length);
        for ( byte b : baseQualities ) {
            convertedBaseQualities.add((int) b);
        }

        genomicsRead.setAlignedQuality(convertedBaseQualities.isEmpty() ? null : convertedBaseQualities);
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

    @Override
    public void setPosition( String contig, int start ) {
        if ( contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1 ) {
            throw new IllegalArgumentException("contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME + ", and start must be >= 1");
        }

        makePositionIfNecessary();

        genomicsRead.getAlignment().getPosition().setReferenceName(contig);
        // Convert from a 1-based to a 0-based position
        genomicsRead.getAlignment().getPosition().setPosition((long) start - 1);
    }

    @Override
    public void setPosition( Locatable locatable ) {
        setPosition(locatable.getContig(), locatable.getStart());
    }

    @Override
    public void setCigar( Cigar cigar ) {
        makeAlignmentIfNecessary();
        genomicsRead.getAlignment().setCigar(CigarConverter.convertSAMCigarToCigarUnitList(cigar));
    }

    @Override
    public void setCigar( String cigarString ) {
        makeAlignmentIfNecessary();
        genomicsRead.getAlignment().setCigar(CigarConverter.convertSAMCigarToCigarUnitList(TextCigarCodec.decode(cigarString)));
    }

    @Override
    public void setIsUnmapped() {
        genomicsRead.setAlignment(null);
    }
}
