package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

/**
 * Class representing the change in the protein sequence for a specific reference/alternate allele pair in a variant.
 * Created by jonn on 10/23/18.
 */
public final class ProteinChangeInfo {

    private final static Logger logger = LogManager.getLogger(ProteinChangeInfo.class);

    /** 1-based inclusive position of the first amino acid changed in this {@link ProteinChangeInfo}. */
    private int    aaStartPos;
    /** 1-based inclusive position of the last amino acid changed in this {@link ProteinChangeInfo}. */
    private int    aaEndPos;
    /**
     * {@link String} representation of the reference amino acid sequence in this {@link ProteinChangeInfo}.
     * May be empty.  Must not be {@code null}.
     */
    private String refAaSeq;
    /**
     * {@link String} representation of the alternate amino acid sequence in this {@link ProteinChangeInfo}.
     * May be empty.  Must not be {@code null}.
     */
    private String altAaSeq;

    private ProteinChangeInfo(final int    aaStartPos,
                              final int    aaEndPos,
                              final String refAaSeq,
                              final String altAaSeq) {
        this.aaStartPos = aaStartPos;
        this.aaEndPos   = aaEndPos;
        this.refAaSeq   = refAaSeq;
        this.altAaSeq   = altAaSeq;
    }

    private ProteinChangeInfo( final Allele refAllele,
                               final Allele altAllele,
                               final int codingSequenceAlleleStart,
                               final int alignedCodingSequenceAlleleStart,
                               final String codingSequence,
                               final Strand strand,
                               final boolean isMitochondria) {

        // Cache whether it's a frameshift variant:
        final boolean isFrameshift =  GATKVariantContextUtils.isFrameshift( refAllele, altAllele );

        // Get our protein sequences:
        final Pair<String, String> proteinSequences = createProteinSequences(refAllele, altAllele, codingSequenceAlleleStart, codingSequence, isFrameshift, isMitochondria);
        final String referenceProteinSequence = proteinSequences.getLeft();
        final String alternateProteinSequence = proteinSequences.getRight();

        // Because we use AminoAcid.UNDECODABLE as a placeholder, we don't actually need to do any more work here.
        // Any UNDECODABLE amino acids will be rendered as their string represntation ("?") in the protein change
        // string.

        // Get the _index_ of the first different amino acid (not the protein position!):
        // Default to the amino acid corresponding to the aligned coding sequence allele start:
        int       proteinChangeStartIndex  = (alignedCodingSequenceAlleleStart - 1) / AminoAcid.CODON_LENGTH;
        final int maxProteinSequenceLength = Math.max(referenceProteinSequence.length(), alternateProteinSequence.length());
        for ( int i = 0; i < maxProteinSequenceLength; ++i ) {
            if ( (i >= referenceProteinSequence.length()) || (i >= alternateProteinSequence.length()) ||
                    (referenceProteinSequence.charAt(i) != alternateProteinSequence.charAt(i)) ) {
                proteinChangeStartIndex = i;
                break;
            }
        }

        // Start pos is correct.
        // Now get the end pos (need to check for FS).

        final boolean indelIsBetweenCodons =
                FuncotatorUtils.isIndelBetweenCodons(
                        codingSequenceAlleleStart,
                        alignedCodingSequenceAlleleStart,
                        refAllele.getBaseString(),
                        strand
                );

        // Get the number of amino acids for which the alt allele codes:
        // Subtract 1 to remove leading base required by VCFs
        final int numAltAminoAcids = (int) Math.ceil((altAllele.length() - 1) / ((double) AminoAcid.CODON_LENGTH));

        // Get the number of amino acids to use as the reference:
        // subtract 1 to remove leading base required by VCFs
        final int numRefAminoAcids = (int) Math.ceil((refAllele.length() - 1) / ((double) AminoAcid.CODON_LENGTH));

        // Frameshifts are always rendered the same way:
        if ( isFrameshift ) {
            initializeForFrameshift(referenceProteinSequence, proteinChangeStartIndex);
        }
        // Handle insertions and deletions:
        else if ( GATKVariantContextUtils.isInsertion(refAllele, altAllele) ) {
            initializeForInsertion(alignedCodingSequenceAlleleStart, strand, referenceProteinSequence, alternateProteinSequence, proteinChangeStartIndex, indelIsBetweenCodons, numAltAminoAcids, numRefAminoAcids);
        }
        else if ( GATKVariantContextUtils.isDeletion(refAllele, altAllele) ) {
            initializeForDeletion(alignedCodingSequenceAlleleStart, strand, referenceProteinSequence, alternateProteinSequence, indelIsBetweenCodons, numAltAminoAcids, numRefAminoAcids);
        }
        else {
            initializeForOnp(referenceProteinSequence, alternateProteinSequence, proteinChangeStartIndex);
        }
    }


    private Pair<String, String> createProteinSequences(final Allele refAllele,
                                                        final Allele altAllele,
                                                        final int codingSequenceAlleleStart,
                                                        final String codingSequence,
                                                        final boolean isFrameshift,
                                                        final boolean isMitochondria) {
        final String referenceProteinSequence;
        final String alternateProteinSequence;

        // Subtract 1 to account for 1-based genomic positions:
        final String altCodingSequence =
                codingSequence.substring(0, codingSequenceAlleleStart - 1) +
                altAllele.getBaseString() +
                codingSequence.substring(codingSequenceAlleleStart + refAllele.length() -1);

        if ( isMitochondria ) {
            // Mitochondrial protein sequences differ from the Standard Code, so we must treat them separately:
            referenceProteinSequence = FuncotatorUtils.createMitochondrialAminoAcidSequence(codingSequence, false, "(size=" + codingSequence.length() + ", ref allele: " + refAllele.getBaseString() + ")");
            alternateProteinSequence = FuncotatorUtils.createMitochondrialAminoAcidSequence( altCodingSequence, isFrameshift, "(size=" + codingSequence.length() + ", alt allele: " + altAllele.getBaseString() + ")");
        }
        else {
            // Create a protein sequence using the Standard Code:
            referenceProteinSequence = FuncotatorUtils.createAminoAcidSequence(codingSequence, false, "(size=" + codingSequence.length() + ", ref allele: " + refAllele.getBaseString() + ")");
            alternateProteinSequence = FuncotatorUtils.createAminoAcidSequence(altCodingSequence, isFrameshift, "(size=" + codingSequence.length() + ", alt allele: " + altAllele.getBaseString() + ")");
        }

        // This part is actually some work to iterate through the string.
        // Only log here if we would warn the user anyway:
        if (logger.isWarnEnabled()) {
            if ( referenceProteinSequence.contains(AminoAcid.UNDECODABLE.getLetter()) ) {
                logger.warn("Ref protein sequence is undecodable: " + codingSequence);
            }
            else if ( alternateProteinSequence.contains(AminoAcid.UNDECODABLE.getLetter()) ) {
                logger.warn("Alt protein sequence is undecodable: " + altCodingSequence);
            }
        }

        return Pair.of(referenceProteinSequence, alternateProteinSequence);
    }

    private void initializeForOnp(final String referenceProteinSequence, final String alternateProteinSequence, final int proteinChangeStartIndex) {
        // ONP - get the length of the change and render the changed bases:
        int i = proteinChangeStartIndex;
        // Go through the protein sequences and find the first Amino Acid that is the same again:
        while (( i < referenceProteinSequence.length()) && (i < alternateProteinSequence.length()) &&
                (referenceProteinSequence.charAt(i) != alternateProteinSequence.charAt(i))) {
            ++i;
        }
        // Get the protein change end index:
        final int proteinChangeEndIndex = i;

        if ( proteinChangeStartIndex == proteinChangeEndIndex ) {
            // We have a single-position / silent change:
            // Add 1 to go from index to 1-based genomic position:
            aaStartPos = proteinChangeStartIndex + 1;
            aaEndPos = aaStartPos;
            refAaSeq = referenceProteinSequence.substring(proteinChangeStartIndex, proteinChangeStartIndex + 1);
            altAaSeq = refAaSeq;
        }
        else {
            // Add 1 to go from index to 1-based genomic position:
            aaStartPos = proteinChangeStartIndex + 1;
            // Leave the end position as we found it to have correct bounds on the Amino Acids:
            aaEndPos = proteinChangeEndIndex;
            refAaSeq = referenceProteinSequence.substring(proteinChangeStartIndex, proteinChangeEndIndex);
            altAaSeq = alternateProteinSequence.substring(proteinChangeStartIndex, proteinChangeEndIndex);
        }
    }

    private void initializeForDeletion(final int alignedCodingSequenceAlleleStart, final Strand strand, final String referenceProteinSequence, final String alternateProteinSequence, final boolean indelIsBetweenCodons, final int numAltAminoAcids, final int numRefAminoAcids) {
        final int proteinChangeStartIndex;// We render the protein change differently if it's a deletion directly between two codons:
        if (indelIsBetweenCodons) {
            // Because we're inbetween codons / have full codons deleted, we can use the start position of the
            // variant:
            // Add 1 to account for the required leading base when on + strands:
            proteinChangeStartIndex = ((alignedCodingSequenceAlleleStart-1) / AminoAcid.CODON_LENGTH) + (strand == Strand.POSITIVE ? 1 : 0);

            aaStartPos = proteinChangeStartIndex + 1;
            aaEndPos = aaStartPos + numRefAminoAcids - 1;

            final int endCoord = Math.min(proteinChangeStartIndex + numRefAminoAcids, referenceProteinSequence.length());

            refAaSeq = referenceProteinSequence.substring(proteinChangeStartIndex, endCoord);
            altAaSeq = "";
        }
        else {
            // To start with, we fill in the information naively corresponding to the potentially
            // changed amino acid sequence:
            proteinChangeStartIndex = ((alignedCodingSequenceAlleleStart - 1) / AminoAcid.CODON_LENGTH);

            // If we're on the - strand, we need to grab 1 fewer amino acid from the end of the sequence:
            final int endOffset = strand == Strand.POSITIVE ? 1 : 0;

            aaStartPos = proteinChangeStartIndex + 1;
            aaEndPos = aaStartPos + numRefAminoAcids + endOffset;

            final int refEndCoord = Math.min(aaEndPos, referenceProteinSequence.length());
            refAaSeq = referenceProteinSequence.substring(proteinChangeStartIndex, refEndCoord);

            final int altEndCoord = Math.min(aaStartPos + numAltAminoAcids + endOffset, alternateProteinSequence.length());
            altAaSeq = alternateProteinSequence.substring(proteinChangeStartIndex, altEndCoord);

            // Trim our state for this deletion:
            trimDeletionProteinChangeVariables();
        }

        // Check to make sure we have any alt amino acids left:
        if ( altAaSeq.isEmpty() ) {
            // We have actually just deleted a set of amino acids.
            // Just set the end to the start to complete the deletion:
            aaEndPos = aaStartPos;
        }
    }

    private void trimDeletionProteinChangeVariables() {
        // OK, now we have a naive, but correct protein change.
        // We should attempt to simplify it
        // into a simple deletion.
        // To do this we must detect the bases that have been deleted from the reference amino acid
        // string.  Because of the properties of in-frame deletions, we know that the middle amino
        // acids are the ones that could be deleted.  Therefore, we can detect this case by checking
        // if first and last amino acids are the same (between refAaSeq and altAaSeq).  If they are
        // then we can safely identify the middle bases as the deleted bases.
        // In the case where only the first or last amino acid is the same, we know that that
        // refAaSeq.length() must be 2, and that it is the other amino acid that has been deleted.
        //
        // Note this case is similar to the complementary insertion case.
        // Note that we must attempte to alternate looking at the front and the back to get the proper
        //      and correct protein change.
        boolean frontMustBeTrimmed = (!refAaSeq.isEmpty()) && (!altAaSeq.isEmpty()) && (refAaSeq.charAt(0) == altAaSeq.charAt(0));
        boolean backMustBeTrimmed  = true;

        while ( frontMustBeTrimmed || backMustBeTrimmed ) {
            if ( frontMustBeTrimmed ) {
                // Set up our data as a simple deletion:
                aaStartPos++;
                aaEndPos = aaStartPos;
                refAaSeq = refAaSeq.substring(1);
                altAaSeq = altAaSeq.substring(1);
            }

            backMustBeTrimmed = (!altAaSeq.isEmpty()) && (refAaSeq.charAt(refAaSeq.length() - 1) == altAaSeq.charAt(altAaSeq.length() - 1));
            if ( backMustBeTrimmed ) {
                // Set up our data as a simple deletion:
                --aaEndPos;
                refAaSeq = refAaSeq.substring(0, refAaSeq.length() - 1);
                altAaSeq = altAaSeq.substring(0, altAaSeq.length() - 1);
            }

            frontMustBeTrimmed = (!refAaSeq.isEmpty()) && (!altAaSeq.isEmpty()) && (refAaSeq.charAt(0) == altAaSeq.charAt(0));
        }
    }

    private void initializeForInsertion(final int alignedCodingSequenceAlleleStart,
                                        final Strand strand,
                                        final String referenceProteinSequence,
                                        final String alternateProteinSequence,
                                        int proteinChangeStartIndex,
                                        final boolean indelIsBetweenCodons,
                                        final int numAltAminoAcids,
                                        final int numRefAminoAcids) {
        // We render the protein change differently if it's an insertion directly between two codons:
        if (indelIsBetweenCodons) {

            // Get the position of the Amino Acid before the insertion:
            aaStartPos = ((alignedCodingSequenceAlleleStart-1) / AminoAcid.CODON_LENGTH) +
                    // If we're on the + strand we need to add 1 to make the amino acid position line up correctly:
                    (strand == Strand.POSITIVE ? 1 : 0);
            aaEndPos = aaStartPos + 1;
            refAaSeq = "";
            final int endCoord = Math.min(proteinChangeStartIndex + numAltAminoAcids, alternateProteinSequence.length());
            altAaSeq = alternateProteinSequence.substring(proteinChangeStartIndex, endCoord );
        }
        else {
            // To start with, we fill in the information naively corresponding to the potentially
            // changed amino acid sequence:
            proteinChangeStartIndex = ((alignedCodingSequenceAlleleStart-1) / AminoAcid.CODON_LENGTH);

            aaStartPos = proteinChangeStartIndex + 1;
            aaEndPos = aaStartPos + numRefAminoAcids;

            refAaSeq = referenceProteinSequence.substring(proteinChangeStartIndex, aaEndPos);
            final int endCoord = Math.min(aaStartPos + numAltAminoAcids, alternateProteinSequence.length());
            altAaSeq = alternateProteinSequence.substring(proteinChangeStartIndex, endCoord);

            // Trim our state for this insertion:
            trimInsertionProteinChangeVariables();
        }
    }

    private void trimInsertionProteinChangeVariables() {
        // OK, now we have a naive, but correct protein change.
        // We should attempt to simplify it into a simple insertion.
        // To do this we must detect the bases that have been inserted into the reference amino acid
        // string.
        // For this case, we check the first and last amino acids in the ref and alt strings to
        // see if they match.  If so we remove that amino acid from the change string.
        //
        // Note this case is similar to the complementary deletion case.
        // Note that we must attempte to alternate looking at the front and the back to get the proper
        //      and correct protein change.
        boolean frontMustBeTrimmed = (!refAaSeq.isEmpty()) && (!altAaSeq.isEmpty()) && (refAaSeq.charAt(0) == altAaSeq.charAt(0));
        boolean backMustBeTrimmed  = true;

        while ( frontMustBeTrimmed || backMustBeTrimmed ) {
            if ( frontMustBeTrimmed ) {
                aaEndPos++;
                refAaSeq = refAaSeq.substring(1);
                altAaSeq = altAaSeq.substring(1);
            }

            backMustBeTrimmed = (!refAaSeq.isEmpty()) && (refAaSeq.charAt(refAaSeq.length() - 1) == altAaSeq.charAt(altAaSeq.length() - 1));
            if ( backMustBeTrimmed ) {

                // Must adjust start position so that insertions can occur with a correct range:
                if ( aaStartPos == aaEndPos ) {
                    --aaStartPos;
                }
                else {
                    --aaEndPos;
                }

                refAaSeq = refAaSeq.substring(0, refAaSeq.length() - 1);
                altAaSeq = altAaSeq.substring(0, altAaSeq.length() - 1);
            }

            frontMustBeTrimmed = (!refAaSeq.isEmpty()) && (!altAaSeq.isEmpty()) && (refAaSeq.charAt(0) == altAaSeq.charAt(0));
        }
    }

    private void initializeForFrameshift(final String referenceProteinSequence, final int proteinChangeStartIndex) {
        aaStartPos = proteinChangeStartIndex + 1;
        aaEndPos = aaStartPos;

        // If we run off the end of the protein without finding a difference, we say the last codon is affected:
        if ( aaEndPos > referenceProteinSequence.length() ) {
            refAaSeq = referenceProteinSequence.substring(referenceProteinSequence.length()-1, referenceProteinSequence.length());
            altAaSeq = "";
        }
        else {
            refAaSeq = referenceProteinSequence.substring(proteinChangeStartIndex, aaEndPos);
            altAaSeq = "";
        }
    }

    /**
     * Create a {@link ProteinChangeInfo} object containing given information about a protein change.
     * @param aaStartPos 1-based inclusive position of the first amino acid changed in this {@link ProteinChangeInfo}.
     * @param aaEndPos 1-based inclusive position of the last amino acid changed in this {@link ProteinChangeInfo}.
     * @param refAaSeq {@link String} representation of the reference amino acid sequence in this {@link ProteinChangeInfo}.  May be empty.  Must not be {@code null}.
     * @param altAaSeq {@link String} representation of the alternate amino acid sequence in this {@link ProteinChangeInfo}.  May be empty.  Must not be {@code null}.
     * @return A new {@link ProteinChangeInfo} object representing the change in the protein sequence for the given input data.
     */
    public static ProteinChangeInfo create(final int    aaStartPos,
                                           final int    aaEndPos,
                                           final String refAaSeq,
                                           final String altAaSeq) {
        return new ProteinChangeInfo(aaStartPos,
                aaEndPos,
                refAaSeq,
                altAaSeq);
    }

    /**
     * Create a {@link ProteinChangeInfo} object which will represent the change in the protein sequence
     * which would be caused by a variant.
     * @param refAllele The strand-corrected (i.e. if on the - strand, it has been reverse-complemented) reference {@link Allele} for the variant.  Must not be {@code null}.
     * @param altAllele The strand-corrected (i.e. if on the - strand, it has been reverse-complemented) alternate {@link Allele} for the variant.  Must not be {@code null}.
     * @param codingSequenceAlleleStart The position (1-based, inclusive) in the _coding sequence_ at which the variant begins.  (NOTE: This is _not_ the same the genomic position, nor is it necessarily the same as the transcript position of the variant).
     * @param alignedCodingSequenceAlleleStart The codon-aligned position (1-based, inclusive) in the _coding sequence_ at which the variant begins.  (NOTE: This is _not_ the same the genomic position, nor is it necessarily the same as the transcript position of the variant).
     * @param codingSequence The strand-corrected (i.e. if on the - strand, it has been reverse-complemented) sequence of bases containing the _coding sequence_ for a particular transcript of a gene, from which we should render a protein change.  (NOTE: This is _not_ the same the gene sequence, nor is it necessarily the same as the whole transcript sequence).  Must not be {@code null}.
     * @param strand The {@link Strand} on which the transcript for this protein change occurs.  Must not be {@link Strand#NONE}.  Must not be {@code null}.
     * @param isMitochondria If {@code true}, will use Mitochondrial protein decoding, rather than the standard eukaryotic amino acid decoding.) {
     * @return A new {@link ProteinChangeInfo} object representing the change in the protein sequence for the given input data.
     */
    public static ProteinChangeInfo create( final Allele refAllele,
                                            final Allele altAllele,
                                            final int codingSequenceAlleleStart,
                                            final int alignedCodingSequenceAlleleStart,
                                            final String codingSequence,
                                            final Strand strand,
                                            final boolean isMitochondria) {
        Utils.nonNull(refAllele);
        Utils.nonNull(altAllele);
        Utils.nonNull(codingSequence);
        Utils.nonNull(strand);

        return new ProteinChangeInfo(refAllele, altAllele, codingSequenceAlleleStart, alignedCodingSequenceAlleleStart, codingSequence, strand, isMitochondria);
    }

    /**
     * @return The amino acid start position (1-based, inclusive) for this {@link ProteinChangeInfo} object.
     */
    public int getAaStartPos() {
        return aaStartPos;
    }

    /**
     * @return The amino acid end position (1-based, inclusive) for this {@link ProteinChangeInfo} object.
     */
    public int getAaEndPos() {
        return aaEndPos;
    }

    /**
     * @return The reference amino acid sequence for this {@link ProteinChangeInfo} object.
     */
    public String getRefAaSeq() {
        return refAaSeq;
    }

    /**
     * @return The alternate amino acid sequence for this {@link ProteinChangeInfo} object.
     */
    public String getAltAaSeq() {
        return altAaSeq;
    }

    @Override
    public boolean equals(final Object that){
        if ( that instanceof ProteinChangeInfo ) {
            final ProteinChangeInfo thatPCI = (ProteinChangeInfo)that;
            return  (aaStartPos == thatPCI.aaStartPos) &&
                    (aaEndPos    == thatPCI.aaEndPos) &&
                    (refAaSeq.equals(thatPCI.refAaSeq)) &&
                    (altAaSeq.equals(thatPCI.altAaSeq));
        }
        return false;
    }

    @Override
    public String toString() {
        return "ProteinChangeInfo{" + aaStartPos + ", " + aaEndPos  + ", " + (refAaSeq.isEmpty() ? "\"\"" : refAaSeq) + ", " + (altAaSeq.isEmpty() ? "\"\"" : altAaSeq) + "}";
    }

    @Override
    public int hashCode() {
        int result = aaStartPos;
        result = 31 * result + aaEndPos;
        result = 31 * result + (refAaSeq != null ? refAaSeq.hashCode() : 0);
        result = 31 * result + (altAaSeq != null ? altAaSeq.hashCode() : 0);
        return result;
    }
}
