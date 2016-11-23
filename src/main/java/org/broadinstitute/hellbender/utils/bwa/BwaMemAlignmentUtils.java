package org.broadinstitute.hellbender.utils.bwa;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Arrays;
import java.util.List;

/**
 * Utils to move data from a BwaMemAlignment into a GATKRead, or into a SAM tag.
 */
public class BwaMemAlignmentUtils {
    /**
     * Transforms an unaligned GATKRead into an aligned one, as specified by the supplied BwaMemAlignment.
     *
     * Also, BWA does this odd thing, supplying AS and XS tags for unaligned reads.
     * To produce SAM records that are byte-identical to the ones created by the bwa mem command line,
     * call this method with alwaysGenerateASandXS=true.
     */
    public static GATKRead applyAlignment(final GATKRead unalignedRead, final BwaMemAlignment alignment,
                                            final List<String> refNames, final SAMFileHeader header,
                                            final boolean softclipAlts, final boolean alwaysGenerateASandXS ) {
        final SAMRecord samRecord = new SAMRecord(header);
        samRecord.setReadName(unalignedRead.getName());
        samRecord.setFlags(alignment.getSamFlag());
        if ( alignment.getRefId() >= 0 ) samRecord.setReferenceName(refNames.get(alignment.getRefId()));
        else if ( alignment.getMateRefId() >= 0 ) samRecord.setReferenceName(refNames.get(alignment.getMateRefId()));
        if ( alignment.getRefStart() >= 0 ) samRecord.setAlignmentStart(alignment.getRefStart()+1);
        else if ( alignment.getMateRefStart() >= 0 ) samRecord.setAlignmentStart(alignment.getMateRefStart()+1);
        if ( alignment.getMapQual() >= 0 ) samRecord.setMappingQuality(alignment.getMapQual());
        byte[] seq = unalignedRead.getBases();
        byte[] quals = unalignedRead.getBaseQualities();
        if ( SAMFlag.READ_REVERSE_STRAND.isSet(alignment.getSamFlag()) &&
                SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(alignment.getSamFlag()) ) {
            seq = BaseUtils.simpleReverseComplement(seq);
            quals = Arrays.copyOf(quals, quals.length);
            SequenceUtil.reverseQualities(quals);
        }
        if ( alignment.getCigar() != null && !alignment.getCigar().isEmpty() ) {
            final Cigar cigar = TextCigarCodec.decode(alignment.getCigar());
            Cigar tmpCigar = cigar;
            if ( !softclipAlts && SAMFlag.SUPPLEMENTARY_ALIGNMENT.isSet(alignment.getSamFlag()) ) {
                if ( tmpCigar.getFirstCigarElement().getOperator() == CigarOperator.S ||
                        tmpCigar.getLastCigarElement().getOperator() == CigarOperator.S ) {
                    tmpCigar = new Cigar();
                    for ( final CigarElement ele : cigar ) {
                        if ( ele.getOperator() == CigarOperator.S ) {
                            tmpCigar.add(new CigarElement(ele.getLength(), CigarOperator.H));
                        }
                        else {
                            tmpCigar.add(ele);
                        }
                    }
                }
                seq = Arrays.copyOfRange(seq, alignment.getSeqStart(), alignment.getSeqEnd());
                quals = Arrays.copyOfRange(quals, alignment.getSeqStart(), alignment.getSeqEnd());
            }
            samRecord.setCigar(tmpCigar);
            samRecord.setAttribute("NM", alignment.getNMismatches());
            samRecord.setAttribute("AS", alignment.getAlignerScore());
            samRecord.setAttribute("XS", alignment.getSuboptimalScore());
            samRecord.setAttribute("MD", alignment.getMDTag());
            samRecord.setAttribute("XA", alignment.getXATag());
        }
        else if ( alwaysGenerateASandXS ) {
            samRecord.setAttribute("AS", 0);
            samRecord.setAttribute("XS", 0);
        }
        if ( alignment.getMateRefId() >= 0 ) samRecord.setMateReferenceName(refNames.get(alignment.getMateRefId()));
        else if ( alignment.getRefId() >= 0 ) samRecord.setMateReferenceName(refNames.get(alignment.getRefId()));
        if ( alignment.getMateRefStart() >= 0 ) samRecord.setMateAlignmentStart(alignment.getMateRefStart()+1);
        else if ( alignment.getRefStart() >= 0 ) samRecord.setMateAlignmentStart(alignment.getRefStart()+1);
        if ( alignment.getTemplateLen() != 0 ) samRecord.setInferredInsertSize(alignment.getTemplateLen());
        else samRecord.setInferredInsertSize(0);
        if ( SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(alignment.getSamFlag()) ) {
            samRecord.setReadBases(seq);
            samRecord.setBaseQualities(quals);
        } else {
            samRecord.setReadBases(SAMRecord.NULL_SEQUENCE);
            samRecord.setBaseQualities(SAMRecord.NULL_QUALS);
        }
        final String readGroup = unalignedRead.getReadGroup();
        if ( readGroup != null ) samRecord.setAttribute(SAMTag.RG.name(), readGroup);
        return SAMRecordToGATKReadAdapter.headerlessReadAdapter(samRecord);
    }

    /**
     * Describes an alignment as a string for use in an SA tag, for example.
     */
    public static String asTag( final BwaMemAlignment alignment, final List<String> refNames ) {
        return refNames.get(alignment.getRefId())+","+(alignment.getRefStart()+1)+","+
                (SAMFlag.READ_REVERSE_STRAND.isSet(alignment.getSamFlag())?"-":"+")+","+
                alignment.getCigar()+","+alignment.getMapQual()+","+alignment.getNMismatches()+";";
    }
}
