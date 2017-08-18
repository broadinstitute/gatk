package org.broadinstitute.hellbender.utils.bwa;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Utils to move data from a BwaMemAlignment into a GATKRead, or into a SAM tag.
 */
public class BwaMemAlignmentUtils {
    /**
     * Builds a SAMRecord from unaligned read data and an alignment.
     * qualsArg can be null.
     * readGroup can be null.
     *
     * Also, BWA does this odd thing, supplying AS and XS tags for unaligned reads.
     * To produce SAM records that are byte-identical to the ones created by the bwa mem command line,
     * call this method with alwaysGenerateASandXS=true.
     */
    public static SAMRecord applyAlignment( final String readName, final byte[] basesArg, final byte[] qualsArg,
                                           final String readGroup, final BwaMemAlignment alignment,
                                           final List<String> refNames, final SAMFileHeader header,
                                           final boolean softclipAlts, final boolean alwaysGenerateASandXS ) {
        Utils.nonNull(alignment, "the input alignment cannot be null");
        Utils.nonNull(basesArg, "the input bases argument cannot be null");
        Utils.nonNull(refNames, "the input reference name list cannot be null");
        final int refId = alignment.getRefId();
        final int refMateId = alignment.getMateRefId();
        if (refId >= refNames.size()) throw new IllegalArgumentException("reference names array to short for refId");
        if (refMateId >= refNames.size()) throw new IllegalArgumentException("reference names array to short for refMateId");
        if (refId >= 0 && refNames.get(refId) == null) throw new IllegalArgumentException("reference name cannot be null");
        if (refMateId >= 0 && refNames.get(refMateId) == null) throw new IllegalArgumentException("mate reference name cannot be null");
        final SAMRecord samRecord = new SAMRecord(header);
        samRecord.setReadName(readName);
        final int samFlag = alignment.getSamFlag();
        samRecord.setFlags(samFlag);
        if ( refId >= 0 ) samRecord.setReferenceName(refNames.get(refId));
        else if ( refMateId >= 0 ) samRecord.setReferenceName(refNames.get(refMateId));
        if ( alignment.getRefStart() >= 0 ) samRecord.setAlignmentStart(alignment.getRefStart()+1);
        else if ( alignment.getMateRefStart() >= 0 ) samRecord.setAlignmentStart(alignment.getMateRefStart()+1);
        if ( alignment.getMapQual() >= 0 ) samRecord.setMappingQuality(alignment.getMapQual());
        byte[] bases = basesArg;
        byte[] quals = qualsArg == null ? new byte[0] : qualsArg;
        if ( SAMFlag.READ_REVERSE_STRAND.isSet(samFlag) &&
                SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(samFlag) ) {
            bases = BaseUtils.simpleReverseComplement(bases);
            quals = Arrays.copyOf(quals, quals.length);
            SequenceUtil.reverseQualities(quals);
        }
        if ( alignment.getCigar() != null && !alignment.getCigar().isEmpty() ) {
            final Cigar cigar = TextCigarCodec.decode(alignment.getCigar());
            Cigar tmpCigar = cigar;
            if ( !softclipAlts && SAMFlag.SUPPLEMENTARY_ALIGNMENT.isSet(samFlag) ) {
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
            }
            if (tmpCigar.getFirstCigarElement().getOperator() == CigarOperator.H ||
                    tmpCigar.getLastCigarElement().getOperator() == CigarOperator.H) {
                final int leftHardClip = tmpCigar.getFirstCigarElement().getOperator() == CigarOperator.H ? tmpCigar.getFirstCigarElement().getLength() : 0;
                final int rightHardClip = tmpCigar.getCigarElements().size() > 1 && tmpCigar.getLastCigarElement().getOperator() == CigarOperator.H ? tmpCigar.getLastCigarElement().getLength() : 0;
                final int firstRightClippedIndex = bases.length - rightHardClip;
                bases = Arrays.copyOfRange(bases, leftHardClip, firstRightClippedIndex);
                if ( quals.length != 0 )
                    quals = Arrays.copyOfRange(quals, leftHardClip, firstRightClippedIndex);
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
        if ( SAMFlag.READ_PAIRED.isSet(samFlag) ) {
            if ( refMateId >= 0 ) samRecord.setMateReferenceName(refNames.get(refMateId));
            else if ( refId >= 0 ) samRecord.setMateReferenceName(refNames.get(refId));
            if ( alignment.getMateRefStart() >= 0 ) samRecord.setMateAlignmentStart(alignment.getMateRefStart() + 1);
            else if ( alignment.getRefStart() >= 0 ) samRecord.setMateAlignmentStart(alignment.getRefStart() + 1);
            if ( alignment.getTemplateLen() != 0 ) samRecord.setInferredInsertSize(alignment.getTemplateLen());
        }
        if ( SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(samFlag) ) {
            samRecord.setReadBases(bases);
            samRecord.setBaseQualities(quals);
        } else {
            samRecord.setReadBases(SAMRecord.NULL_SEQUENCE);
            samRecord.setBaseQualities(SAMRecord.NULL_QUALS);
        }
        if ( readGroup != null ) samRecord.setAttribute(SAMTag.RG.name(), readGroup);
        return samRecord;
    }

    /**
     * Produces an SA tag for each primary line and supplemental alignment as a map from alignment to tag value.
     */
    public static Map<BwaMemAlignment,String> createSATags( final List<BwaMemAlignment> alignments,
                                                            final List<String> refNames ) {
        final int nAlignments = alignments.size();
        if ( nAlignments < 2 ) return Collections.emptyMap();

        final String[] selfTags = new String[nAlignments];
        for ( int idx = 0; idx != nAlignments; ++idx ) {
            final BwaMemAlignment alignment = alignments.get(idx);
            if ( SAMFlag.NOT_PRIMARY_ALIGNMENT.isUnset(alignment.getSamFlag()) ) {
                selfTags[idx] = asTag(alignment, refNames);
            }
        }

        final Map<BwaMemAlignment,String> saTags = new HashMap<>(SVUtils.hashMapCapacity(nAlignments));
        for ( int idx = 0; idx != nAlignments; ++idx ) {
            final BwaMemAlignment alignment = alignments.get(idx);
            if ( SAMFlag.NOT_PRIMARY_ALIGNMENT.isSet(alignment.getSamFlag()) ) continue;

            final StringBuilder saTag = new StringBuilder();
            for ( int idx2 = 0; idx2 != nAlignments; ++idx2 ) {
                if ( idx2 == idx ) continue;
                final String tag = selfTags[idx2];
                if ( tag != null ) saTag.append(tag);
            }
            saTags.put(alignment, saTag.toString());
        }
        return saTags;
    }

    /**
     * Describes an alignment as a string for use in an SA tag, for example.
     */
    public static String asTag( final BwaMemAlignment alignment, final List<String> refNames ) {
        return refNames.get(alignment.getRefId())+","+(alignment.getRefStart()+1)+","+
                (SAMFlag.READ_REVERSE_STRAND.isSet(alignment.getSamFlag())?"-":"+")+","+
                alignment.getCigar().replace('H','S')+","+alignment.getMapQual()+","+alignment.getNMismatches()+";";
    }
}
