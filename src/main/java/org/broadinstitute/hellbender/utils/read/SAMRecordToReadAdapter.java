package org.broadinstitute.hellbender.utils.read;


import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;

public final class SAMRecordToReadAdapter implements MutableRead {

    private final SAMRecord samRecord;

    public SAMRecordToReadAdapter( final SAMRecord samRecord ) {
        this.samRecord = samRecord;
    }

    @Override
    public String getContig() {
        if ( isUnmapped() ) {
            return SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
        }

        return samRecord.getReferenceName();
    }

    @Override
    public int getStart() {
        if ( isUnmapped() ) {
            return GenomicIndexUtil.UNSET_GENOMIC_LOCATION;
        }

        return samRecord.getAlignmentStart();
    }

    @Override
    public int getEnd() {
        if ( isUnmapped() ) {
            return GenomicIndexUtil.UNSET_GENOMIC_LOCATION;
        }

        return samRecord.getAlignmentEnd();
    }

    @Override
    public String getName() {
        return samRecord.getReadName();
    }

    @Override
    public int getLength() {
        return samRecord.getReadLength();
    }

    @Override
    public byte[] getBases() {
        final byte[] bases = samRecord.getReadBases();
        return bases != null ? bases : new byte[0];
    }

    @Override
    public byte[] getBaseQualities() {
        final byte[] baseQualities = samRecord.getBaseQualities();
        return baseQualities != null ? baseQualities : new byte[0];
    }

    @Override
    public Cigar getCigar() {
        final Cigar cigar = samRecord.getCigar();
        return cigar != null ? cigar : new Cigar();
    }

    @Override
    public boolean isUnmapped() {
        return samRecord.getReadUnmappedFlag() ||
               samRecord.getReferenceIndex() == null || samRecord.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX ||
               samRecord.getReferenceName() == null || samRecord.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
               samRecord.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START;
    }

    @Override
    public void setName( String name ) {
        samRecord.setReadName(name);
    }

    @Override
    public void setBases( byte[] bases ) {
        samRecord.setReadBases(bases);
    }

    @Override
    public void setBaseQualities( byte[] baseQualities ) {
        samRecord.setBaseQualities(baseQualities);
    }

    @Override
    public void setPosition( String contig, int start ) {
        // Don't allow client to set read to be unmapped via this method
        if ( contig == null || contig.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) || start < 1 ) {
            throw new IllegalArgumentException("contig must be non-null and not equal to " + SAMRecord.NO_ALIGNMENT_REFERENCE_NAME + ", and start must be >= 1");
        }

        samRecord.setReferenceName(contig);
        samRecord.setAlignmentStart(start);
        samRecord.setReadUnmappedFlag(false);
    }

    @Override
    public void setPosition( Locatable locatable ) {
        setPosition(locatable.getContig(), locatable.getStart());
    }

    @Override
    public void setCigar( Cigar cigar ) {
        samRecord.setCigar(cigar);
    }

    @Override
    public void setCigar( String cigarString ) {
        samRecord.setCigarString(cigarString);
    }

    @Override
    public void setIsUnmapped() {
        samRecord.setReadUnmappedFlag(true);
    }
}