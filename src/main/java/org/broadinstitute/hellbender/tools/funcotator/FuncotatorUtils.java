package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

public final class FuncotatorUtils {

    private final static Logger logger = Logger.getLogger(FuncotatorUtils.class);

    /**
     * PRIVATE CONSTRUCTOR
     * DO NOT INSTANTIATE THIS CLASS!
     */
    private FuncotatorUtils() {}

    private static final Map<String, AminoAcid> tableByCodon;
    private static final Map<String, AminoAcid> tableByCode;
    private static final Map<String, AminoAcid> tableByLetter;

    /**
     * Initialize our hashmaps of lookup tables:
     */
    static {

        final HashMap<String, AminoAcid> mapByCodon  = new HashMap<>(AminoAcid.values().length);
        final HashMap<String, AminoAcid> mapByLetter = new HashMap<>(AminoAcid.values().length);
        final HashMap<String, AminoAcid> mapByCode   = new HashMap<>(AminoAcid.values().length);

        for ( final AminoAcid acid : AminoAcid.values() ) {
            mapByCode.put(acid.getCode(),acid);
            mapByLetter.put(acid.getLetter(), acid);
            for ( final String codon : acid.getCodons() ) {
                mapByCodon.put(codon,acid);
            }
        }

        tableByCodon = Collections.unmodifiableMap(mapByCodon);
        tableByCode = Collections.unmodifiableMap(mapByCode);
        tableByLetter = Collections.unmodifiableMap(mapByLetter);
    }

    /**
     * Returns the {@link AminoAcid} corresponding to the given three-letter Eukaryotic {@code codon}
     * The codons given are expected to be valid for Eukaryotic DNA.
     * @param codon The three-letter codon (each letter one of A,T,G,C) representing a Eukaryotic {@link AminoAcid}
     * @return The {@link AminoAcid} corresponding to the given {@code codon}.  Returns {@code null} if the given {@code codon} does not code for a Eucaryotic {@link AminoAcid}.
     */
    public static AminoAcid getEukaryoticAminoAcidByCodon(final String codon) {
        if (codon == null) {
            return null;
        }
        return tableByCodon.get(codon.toUpperCase());
    }

    /**
     * Returns the {@link AminoAcid} corresponding to the given three-letter Mitochondrial {@code codon}.
     * The codons given are expected to be valid for Mitochondrial DNA.
     * @param codon The three-letter codon (each letter one of A,T,G,C) representing a Mitochondrial {@link AminoAcid}
     * @return The {@link AminoAcid} corresponding to the given {@code codon}.  Returns {@code null} if the given {@code codon} does not code for a Mitochondrial {@link AminoAcid}.
     */
    public static AminoAcid getMitochondrialAminoAcidByCodon(final String codon, final boolean isFirst) {

        if (codon == null) {
            return null;
        }

        final String upperCodon = codon.toUpperCase();
        if ( isFirst && upperCodon.equals("ATT") || upperCodon.equals("ATA") ) {
            return AminoAcid.METHIONINE;
        } else if ( upperCodon.equals("AGA") || upperCodon.equals("AGG") ) {
            return AminoAcid.STOP_CODON;
        } else if ( upperCodon.equals("TGA") ) {
            return AminoAcid.TRYPTOPHAN;
        } else {
            return tableByCodon.get(upperCodon);
        }
    }

    /**
     * Returns the {@link AminoAcid} corresponding to the given single-letter abbreviation.
     * @param letter The one-letter abbreviation representing an {@link AminoAcid}
     * @return The {@link AminoAcid} corresponding to the given {@code letter}.  Returns {@code null} if the given {@code letter} does not code for an {@link AminoAcid}.
     */
    public static AminoAcid getAminoAcidByLetter(final String letter) {
        if ( letter == null ) {
            return null;
        }

        return tableByLetter.get(letter);
    }

    /**
     * Returns the {@link AminoAcid} corresponding to the given single-letter abbreviation.
     * @param letter The one-letter abbreviation representing an {@link AminoAcid}
     * @return The {@link AminoAcid} corresponding to the given {@code letter}.  Returns {@code null} if the given {@code letter} does not code for an {@link AminoAcid}.
     */
    public static AminoAcid getAminoAcidByLetter(final char letter) {
        return tableByLetter.get(String.valueOf(letter));
    }

    /**
     * @return A {@link String} array of long names for all amino acids in {@link AminoAcid}
     */
    public static String[] getAminoAcidNames() {
        final String[] names = new String[AminoAcid.values().length];
        for ( final AminoAcid acid : AminoAcid.values() ) {
            names[acid.ordinal()] = acid.getName();
        }

        return names;
    }

    /**
     * @return A {@link String} array of short names / three-letter abbreviations for all amino acids in {@link AminoAcid}
     */
    public static String[] getAminoAcidCodes() {
        final String[] codes = new String[AminoAcid.values().length];
        for ( final AminoAcid acid : AminoAcid.values() ) {
            codes[acid.ordinal()] = acid.getCode();
        }

        return codes;
    }

    /**
     * Gets whether the given position is in-frame.
     * @param position The position (1-based, inclusive) to check for alignment / in-frameness.  Assumed to be in a coding-region.  Must be > 0;
     * @return {@code true} if the given {@code position} is in-frame; {@code false} otherwise.
     */
    public static boolean isPositionInFrame( final int position ) {

        ParamUtils.isPositive(position, "Genomic positions start at 1.");

        return (((position - 1) % 3) == 0);
    }

    /**
     * Get the string of bases that are different from the given alleles.
     * Assumes there is one contiguous string of changed bases between the two alleles.
     * Assumes that if there is overlap between the alleles, the overlap occurs at either the front or the back.
     * @param firstAllele First {@link Allele}.  Must not be {@code null}.
     * @param secondAllele Second {@link Allele}.  Must not be {@code null}.
     * @param copyRefBasesWhenAltIsPastEnd Will copy the bases from the given {@code firstAllele} when the alternate allele has no more bases to copy over.  Used primarily for handling deletions.
     * @return A string containing the bases from the given {@code secondAllele} that are different from the reference (in their correct relative order).
     */
    public static String getNonOverlappingAltAlleleBaseString( final Allele firstAllele, final Allele secondAllele, final boolean copyRefBasesWhenAltIsPastEnd ) {

        final StringBuilder sb = new StringBuilder();

        final char noBase = 'x';
        final int maxAlleleLength = Math.max(firstAllele.length(), secondAllele.length());
        boolean havePassedOverlapAlready = false;

        // Find out where the overlap is:
        if ( firstAllele.getBases()[0] != secondAllele.getBases()[0] ) {
            // overlap in back:
            for ( int i = 0; i < maxAlleleLength; ++i ) {
                // Check for differences between the bases:
                char firstAlleleBase = noBase;
                char secondAlleleBase = noBase;

                if ( (firstAllele.length() - 1 - i) >= 0 ) {
                    firstAlleleBase = firstAllele.getBaseString().charAt(firstAllele.length() - 1 - i);
                }
                if ( (secondAllele.length() - 1 - i) >= 0 ) {
                    secondAlleleBase = secondAllele.getBaseString().charAt(secondAllele.length() - 1 - i);
                }

                // Check to see if we're at the end of the differences between the alleles:
                if ( secondAlleleBase == noBase ) {
                    if ( copyRefBasesWhenAltIsPastEnd ) {
                        sb.append( firstAllele.getBaseString().substring(0, firstAllele.length() - i) );
                    }
                    break;
                }
                else if ( havePassedOverlapAlready && (secondAlleleBase == firstAlleleBase) ) {
                    // We can now safely copy the rest of the alt allele into the string buffer and get out of the loop.
                    sb.append( secondAllele.getBaseString().substring(0, secondAllele.length() - i) );
                    break;
                }
                else if (secondAlleleBase != firstAlleleBase) {
                    sb.append(secondAlleleBase);
                    havePassedOverlapAlready = true;
                }
            }

            // Now we rotate the string buffer because we've been iterating through it backwards:
            sb.reverse();
        }
        else {
            // overlap in front:
            for ( int i = 0; i < maxAlleleLength; ++i ) {

                // Check for differences between the bases:
                char firstAlleleBase = noBase;
                char secondAlleleBase = noBase;

                if ( i < firstAllele.length() ) {
                    firstAlleleBase = firstAllele.getBaseString().charAt(i);
                }
                if ( i < secondAllele.length() ) {
                    secondAlleleBase = secondAllele.getBaseString().charAt(i);
                }

                // Check to see if we're at the end of the differences between the alleles:
                if ( secondAlleleBase == noBase ) {
                    if ( copyRefBasesWhenAltIsPastEnd ) {
                        sb.append( firstAllele.getBaseString().substring(i) );
                    }
                    break;
                }
                else if ( havePassedOverlapAlready && (secondAlleleBase == firstAlleleBase) ) {
                    // We can now safely copy the rest of the alt allele into the string buffer and get out of the loop.
                    sb.append( secondAllele.getBaseString().substring(i) );
                    break;
                }
                else if ( secondAlleleBase != firstAlleleBase ) {
                    sb.append(secondAlleleBase);
                    havePassedOverlapAlready = true;
                }
            }
        }

        // We're done.  Bye bye!
        return sb.toString();
    }

    /**
     * Gets the position describing where the given allele and variant lie inside the given transcript using transcript-based coordinates.
     * The index will be calculated even if the given variant ends outside the bounds of the given transcript.
     * Assumes {@code transcript} is a sorted list (in exon number order).
     * @param variant A {@link Locatable} to locate inside the given {@code transcript}.  Must not be {@code null}.
     * @param transcript A sorted {@link List} of {@link Locatable} (in exon number order) that describe the transcript to use for locating the given {@code allele}.  Must be on the same {@code contig} as the given {@code variant}.  Must not be {@code null}.
     * @param strand The strand on which the transcript is read.    Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return The position (1-based, inclusive) describing where the given {@code allele} lies in the given {@code transcript}.  If the variant is not in the given {@code transcript}, then this returns -1.
     */
    public static int getStartPositionInTranscript( final Locatable variant,
                                                    final List<? extends Locatable> transcript,
                                                    final Strand strand) {
        Utils.nonNull(variant);
        Utils.nonNull(transcript);
        assertValidStrand( strand );

        int position = 1;

        boolean foundPosition = false;

        // Creating a new SimpleInterval here so that the code can have fewer checks for strandedness:
        final SimpleInterval variantStartLocus;
        if ( strand == Strand.POSITIVE ) {
            variantStartLocus = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getStart());
        }
        else {
            variantStartLocus = new SimpleInterval(variant.getContig(), variant.getEnd(), variant.getEnd());
        }

        // Ensure all exons are on the same contig as the variant:
        if ( !transcript.stream().allMatch(exon -> exon.getContig().equals(variantStartLocus.getContig())) ) {
            throw new GATKException("Variant and transcript contigs are not equal.");
        }

        for (final Locatable exon : transcript) {

            if (new SimpleInterval(exon).contains(variantStartLocus)) {
                if ( strand == Strand.POSITIVE ) {
                    position += variantStartLocus.getStart() - exon.getStart();
                }
                else {
                    position += exon.getEnd() - variantStartLocus.getStart();
                }
                foundPosition = true;
                break;
            } else {
                // Add 1 because of inclusive positions / indexing starting at 1
                position += exon.getEnd() - exon.getStart() + 1;
            }
        }

        if ( foundPosition ) {
            return position;
        }

        return -1;
    }

    /**
     * Get the sequence-aligned end position for the given allele end position.
     * @param alleleEndPosition The genome end position (1-based, inclusive) for an allele.  Must be > 0.
     * @return An aligned end position (1-based, inclusive) for the given allele end position.
     */
    public static int getAlignedEndPosition(final int alleleEndPosition) {

        ParamUtils.isPositive( alleleEndPosition, "Genomic positions must be > 0." );

        return (int)(Math.ceil(alleleEndPosition / 3.0) * 3);
    }

    /**
     * Gets the sequence aligned position (1-based, inclusive) for the given coding sequence position.
     * This will produce the next lowest position evenly divisible by 3, such that a codon starting at this returned
     * position would include the given position.
     * @param position A sequence starting coordinate for which to produce an coding-aligned position.  Must be > 0.
     * @return A coding-aligned position (1-based, inclusive) corresponding to the given {@code position}.
     */
    public static int getAlignedPosition(final int position) {

        ParamUtils.isPositive( position, "Genomic positions must be > 0." );

        return position - ((position - 1) % 3);
    }

    /**
     * Calculates whether the given {@code startPosition} (1-based, inclusive) is in frame relative to the end of the region.
     * @param startPosition The position (1-based, inclusive) relative to the start of a region to check for frame alignment.    Must be > 0.
     * @param regionLength The length of the region containing {@code startPosition}.  Must be >= 0.
     * @return {@code true} if the given {@code startPosition} is in frame relative to the given {@code regionLength} ; {@code false} otherwise.
     */
    public static boolean isInFrameWithEndOfRegion(final int startPosition, final int regionLength) {

        ParamUtils.isPositive( startPosition, "Genomic positions must be > 0." );
        ParamUtils.isPositiveOrZero( regionLength, "Region length must be >= 0." );

        return (((regionLength - startPosition + 1) % 3) == 0);
    }

    /**
     * Checks to see whether a given indel location occurs on a codon boundary.
     * That is, whether the given indel location is not within a codon but is cleanly between two adjacent codons.
     * NOTE: ASSUMES that there is a leading base that is not part of the indel prepended to the indel string for context.
     * @param codingSequenceAlleleStart The start position of the variant in the coding sequence.
     * @param alignedCodingSequenceAlleleStart The start position of the first codon containing part of the variant in the coding sequence.
     * @param refAllele A {@link String} containing the bases in the reference allele for the given variant.
     * @return {@code true} if the given indel cleanly occurs between two adjacent codons; {@code false} otherwise.
     */
    private static boolean isIndelBetweenCodons(final int codingSequenceAlleleStart,
                                                final int alignedCodingSequenceAlleleStart,
                                                final String refAllele ) {

        final int codonOffset = codingSequenceAlleleStart - alignedCodingSequenceAlleleStart;
        return (((codonOffset + refAllele.length()) % 3) == 0);
    }

    /**
     * Create a properly capitalized version of the given alternate allele codon change string for an insertion by
     * capitalizing only the bases that are different from the reference.
     * This is done using the index of where the alternate allele starts and the length of the alternate allele.
     * @param alternateAlleleCodonChangeString The alternate allele codon {@link String} to capitalize properly.
     * @param alternateAllele A {@link String} containing the bases in the alternate allele.
     * @param startingOffset The offset in {@code alternateAlleleCodonChangeString} where the alternate allele begins.
     * @return A properly capitalized version of the given {@code alternateAlleleCodonChangeString} based on the given alternate allele information.
     */
    private static String createCapitalizedAlternateAlleleInsertionCodonChangeString(final String alternateAlleleCodonChangeString,
                                                                                     final String alternateAllele,
                                                                                     final int startingOffset ) {
        final StringBuilder sb = new StringBuilder();

        // Account for the leading base that we require for insertions:
        final int newStartingOffset = startingOffset + 1;

        int i = 0;
        while ( i < newStartingOffset ) {
            sb.append(Character.toLowerCase( alternateAlleleCodonChangeString.charAt(i++) ) );
        }
        // Subtract 1 from the allele length to account for the leading base that we require for insertions:
        while ( (i - newStartingOffset) < (alternateAllele.length() - 1) ) {
            sb.append(Character.toUpperCase( alternateAlleleCodonChangeString.charAt(i++) ) );
        }
        while ( i < alternateAlleleCodonChangeString.length() ) {
            sb.append(Character.toLowerCase( alternateAlleleCodonChangeString.charAt(i++) ) );
        }

        return sb.toString();
    }

    /**
     * Creates the string representation of the codon change for the given {@link SequenceComparison}.
     * Requires that the given {@code seqComp} has the following fields defined with values that are not {@code null}:
     *     alignedCodingSequenceAlleleStart
     *     alignedReferenceAlleleStop
     *     referenceAllele
     *     alignedCodingSequenceReferenceAllele
     *     alternateAllele
     *     alignedCodingSequenceAlternateAllele
     *     codingSequenceAlleleStart
     * @param seqComp {@link SequenceComparison} representing the alternate and reference alleles for a DNA sequence.  Must not be {@code null}.
     * @return A {@link String} representing the codon change for the given {@link SequenceComparison}.
     */
    public static String getCodonChangeString( final SequenceComparison seqComp ) {

        // ONP:
        if ( GATKVariantContextUtils.isXnp(seqComp.getAlignedReferenceAllele(), seqComp.getAlignedAlternateAllele()) ) {
            return getCodonChangeStringForOnp(
                    seqComp.getAlignedCodingSequenceReferenceAllele(),
                    seqComp.getAlignedCodingSequenceAlternateAllele(),
                    seqComp.getAlignedCodingSequenceAlleleStart(),
                    seqComp.getAlignedReferenceAlleleStop()
            );
        }
        else {

            // Insertion:
            if ( GATKVariantContextUtils.isInsertion(seqComp.getAlignedReferenceAllele(), seqComp.getAlignedAlternateAllele()) ) {

                if ( GATKVariantContextUtils.isFrameshift(seqComp.getAlignedReferenceAllele(), seqComp.getAlignedAlternateAllele()) ) {

                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {
                        final String nextRefCodon = getNextReferenceCodon(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());
                        return "c.(" + (seqComp.getAlignedCodingSequenceAlleleStart() + 3) + "-" + (seqComp.getAlignedCodingSequenceAlleleStart() + 5) + ")" +
                                nextRefCodon.toLowerCase() + "fs";
                    }
                    else {

                        // Right now, we don't care about anything that comes after the frameshift because it'll be garbage.
                        // Therefore, since we know that we're "inside" the effected codon, we report it as the frameshifted
                        // codon:
                        return "c.(" + (seqComp.getAlignedCodingSequenceAlleleStart()) + "-" + (seqComp.getAlignedCodingSequenceAlleleStart() + 2) + ")" +
                                seqComp.getAlignedCodingSequenceReferenceAllele().toLowerCase() + "fs";
                    }
                }
                else {
                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {
                        final String nextRefCodon = getNextReferenceCodon(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());

                        // Here we adjust everything "right" by 3 bases (1 codon) because of the leading base that is
                        // required for indels:
                        return "c.(" + (seqComp.getAlignedCodingSequenceAlleleStart() + 3) + "-" + (seqComp.getAlignedCodingSequenceAlleleStart() + 5) + ")" +
                                nextRefCodon.toLowerCase() + ">" + seqComp.getAlignedAlternateAllele().substring(3) + nextRefCodon.toLowerCase();
                    }
                    else {

                        // Get the number of extra codons to get:
                        // NOTE: We only ever get 1 extra codon if we might overlap with the next codon
                        //       (as per oncotator conventions):
                        final int numAdditionalCodonsToGet = ((seqComp.getCodingSequenceAlleleStart() + 1) % 3 == 0 ? 1 : 0);

                        // Get the next few required reference codons:
                        final List<String> nextRefCodons = getNextReferenceCodons(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand(), numAdditionalCodonsToGet);

                        final String capitalizedAltString = createCapitalizedAlternateAlleleInsertionCodonChangeString(
                                seqComp.getAlignedAlternateAllele() + String.join("", nextRefCodons).toLowerCase(),
                                seqComp.getAlternateAllele(),
                                seqComp.getCodingSequenceAlleleStart() - seqComp.getAlignedCodingSequenceAlleleStart()
                        );

                        return "c.(" + seqComp.getAlignedCodingSequenceAlleleStart() + "-" + (seqComp.getAlignedCodingSequenceAlleleStart() + 2 + (numAdditionalCodonsToGet * 3)) + ")" +
                                seqComp.getAlignedCodingSequenceReferenceAllele().toLowerCase() + String.join("", nextRefCodons).toLowerCase() + ">" + capitalizedAltString;
                    }
                }

            }
            // Deletion:
            else {
                if ( GATKVariantContextUtils.isFrameshift(seqComp.getAlignedReferenceAllele(), seqComp.getAlignedAlternateAllele()) ) {
                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {

                        // Check to see if the deletion actually starts in the next codon, if it does then we skip the
                        // current codon and only display the next one:
                        if ((seqComp.getCodingSequenceAlleleStart() % 3) == 0) {
                            final String nextRefCodon = getNextReferenceCodon(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());
                            return "c.(" + (seqComp.getAlignedCodingSequenceAlleleStart() + 3) + "-" + (seqComp.getAlignedCodingSequenceAlleleStart() + 5) + ")" +
                                    nextRefCodon.toLowerCase() + "fs";
                        }
                        else {
                            return "c.(" + seqComp.getAlignedCodingSequenceAlleleStart() + "-" + seqComp.getAlignedReferenceAlleleStop() + ")" +
                                    seqComp.getAlignedCodingSequenceReferenceAllele().toLowerCase() + "fs";
                        }
                    }
                    else {

                        // We started in the middle of a codon.
                        // We should report the current codon as the start position.

                        // Check to see if the deletion actually starts in the next codon, if it does then we skip the
                        // current codon and only display the next one:
                        if ((seqComp.getCodingSequenceAlleleStart() % 3) == 0) {
                            // First base (i.e. the required base before indels) is the last base of the current
                            // codon.  So we can skip it.
                            return "c.(" + (seqComp.getAlignedCodingSequenceAlleleStart() + 3) + "-" + seqComp.getAlignedReferenceAlleleStop() + ")" +
                                    seqComp.getAlignedCodingSequenceReferenceAllele().substring(3).toLowerCase() + "fs";
                        }
                        else {

                            return "c.(" + (seqComp.getAlignedCodingSequenceAlleleStart()) + "-" + seqComp.getAlignedReferenceAlleleStop() + ")" +
                                    seqComp.getAlignedCodingSequenceReferenceAllele().toLowerCase() + "fs";
                        }
                    }
                }
                else {

                    // Determine how many codons to get:
                    final int numAdditionalCodonsToGet = (int) Math.ceil((seqComp.getAlternateAllele().length() + ((seqComp.getAlignedCodingSequenceAlleleStart() - seqComp.getCodingSequenceAlleleStart()) % 3)) / 3);

                    // Get the next few required reference codons:
                    final List<String> nextRefCodons = getNextReferenceCodons(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand(), numAdditionalCodonsToGet);

                    // Get the end position of the variant:
                    final int endPos = seqComp.getAlignedReferenceAlleleStop() + (3 * numAdditionalCodonsToGet);

                    // Get the string for all reference codons (including the next few codons we need):
                    final String allRefCodons = (seqComp.getAlignedCodingSequenceReferenceAllele() + String.join("", nextRefCodons)).toLowerCase();

                    // This means that the deletion is aligned with a codon:
                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {

                        // An entire codon / set of codons got deleted.
                        // Just put `DEL` after it.
                        // NOTE: We need to ignore the first codon in the list because of the leading base for indels
                        return "c.(" + (seqComp.getAlignedCodingSequenceAlleleStart() + 3) + "-" + endPos
                                + ")" + allRefCodons.substring(3) + "del";
                    }
                    // Non-frameshift deletion starting within a codon:
                    else {
                        // We must report all old codons and then the new codons.
                        // However, for the new codons, we need to create a codon string that includes the new alternate
                        // allele and ONLY the bases that follow it (rounding to the next codon boundary).
                        return "c.(" + seqComp.getAlignedCodingSequenceAlleleStart() + "-" + endPos
                                + ")" + allRefCodons + ">" + seqComp.getAlignedCodingSequenceAlternateAllele().toLowerCase() + String.join("", nextRefCodons).toLowerCase();
                    }
                }
            }
        }
    }

    /**
     * Get the codon change string for an ONP.
     * @param alignedRefAllele The {@link String} of bases contained in the Reference allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param alignedAltAllele The {@link String} of bases contained in the Alternate allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param alignedCodingSequenceAlleleStart The position (1-based, inclusive) of the start of the of the first of the alleles in the coding sequence (prepended with up to 2 bases to be in-frame).
     * @param alignedReferenceAlleleStop The position (1-based, inclusive) of the in-frame end of the Reference allele.
     * @return The {@link String} containing the codon change for the given alleles.
     */
    private static String getCodonChangeStringForOnp( final String alignedRefAllele,
                                                      final String alignedAltAllele,
                                                      final int alignedCodingSequenceAlleleStart,
                                                      final int alignedReferenceAlleleStop){
        final StringBuilder ref = new StringBuilder();
        final StringBuilder alt = new StringBuilder();

        // Capitalize the right parts of each string if they're of equal length:
        for (int i = 0; i < alignedRefAllele.length(); ++i) {
            if (alignedRefAllele.charAt(i) != alignedAltAllele.charAt(i)) {
                ref.append(Character.toUpperCase(alignedRefAllele.charAt(i)));
                alt.append(Character.toUpperCase(alignedAltAllele.charAt(i)));
            } else {
                final char c = Character.toLowerCase(alignedRefAllele.charAt(i));
                ref.append(c);
                alt.append(c);
            }
        }

        // Construct and return the string:
        if (alignedCodingSequenceAlleleStart == alignedReferenceAlleleStop) {
            return "c.(" + alignedCodingSequenceAlleleStart + ")" +
                    ref.toString() + ">" + alt.toString();
        } else {
            return "c.(" + alignedCodingSequenceAlleleStart + "-" +
                    alignedReferenceAlleleStop + ")" +
                    ref.toString() + ">" + alt.toString();
        }
    }

    /**
     * Gets the next complete in-frame codon from the given {@link ReferenceSequence} according to the current codon position and strand.
     * @param referenceSequence The {@link ReferenceSequence} for the current codon.  Must not be {@code null}.
     * @param currentAlignedCodingSequenceAlleleStart The starting position (1-based, inclusive) of the current codon.  Must be > 0.
     * @param currentAlignedCodingSequenceAlleleStop The ending position (1-based, inclusive) of the current codon.  Must be > 0.
     * @param strand The {@link Strand} on which the current codon resides.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return The next codon in frame with the current codon as specified by the given current codon positions.
     */
    private static String getNextReferenceCodon(final ReferenceSequence referenceSequence,
                                                final int currentAlignedCodingSequenceAlleleStart,
                                                final int currentAlignedCodingSequenceAlleleStop,
                                                final Strand strand) {

        Utils.nonNull( referenceSequence );
        ParamUtils.isPositive(currentAlignedCodingSequenceAlleleStart, "Genomic positions must be > 0.");
        ParamUtils.isPositive(currentAlignedCodingSequenceAlleleStop, "Genomic positions must be > 0.");
        assertValidStrand(strand);

        final String nextRefCodon;
        if ( strand == Strand.POSITIVE ) {
            nextRefCodon = referenceSequence.getBaseString().substring(currentAlignedCodingSequenceAlleleStop, currentAlignedCodingSequenceAlleleStop + 3 );
        }
        else {
            nextRefCodon = ReadUtils.getBasesReverseComplement(
                    referenceSequence.getBaseString().substring(currentAlignedCodingSequenceAlleleStart - 3, currentAlignedCodingSequenceAlleleStart ).getBytes()
            );
        }
        return nextRefCodon;
    }

    /**
     * Gets the requested number of complete in-frame codons from the given {@link ReferenceSequence} that follow the given current codon position and strand.
     * @param referenceSequence The {@link ReferenceSequence} for the current codon.  Must not be {@code null}.
     * @param currentAlignedCodingSequenceAlleleStart The starting position (1-based, inclusive) of the current codon.  Must be > 0.
     * @param currentAlignedCodingSequenceAlleleStop The ending position (1-based, inclusive) of the current codon.  Must be > 0.
     * @param strand The {@link Strand} on which the current codon resides.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @param numAdditionalCodonsToGet The number of codons to return.
     * @return The {@link List} of codons (as {@link String}s) in frame with the current codon as specified by the given current codon positions.
     */
    private static List<String> getNextReferenceCodons(final ReferenceSequence referenceSequence,
                                                       final int currentAlignedCodingSequenceAlleleStart,
                                                       final int currentAlignedCodingSequenceAlleleStop,
                                                       final Strand strand,
                                                       final int numAdditionalCodonsToGet) {

        ParamUtils.isPositiveOrZero( numAdditionalCodonsToGet, "Must specify a positive number of codons to return (or zero)." );

        final ArrayList<String> nextCodons = new ArrayList<>(numAdditionalCodonsToGet);

        for (int i = 0; i < numAdditionalCodonsToGet; ++i) {
            final String nextCodon = getNextReferenceCodon(referenceSequence, currentAlignedCodingSequenceAlleleStart + (i*3), currentAlignedCodingSequenceAlleleStop + (i*3), strand);
            nextCodons.add(nextCodon);
        }

        return nextCodons;
    }

    /**
     * Gets a codon change string for a splice site.
     * Assumes the variant and exon referenced in the params are on the same contig.
     * @param variantStart Start position (1-based, inclusive) of the variant.  Must be > 0.
     * @param exonNumber Number (1-based) of the exon in the transcript.  Must be > 0.
     * @param exonStart Start position (1-based, inclusive) of the exon.  Must be > 0.
     * @param exonEnd End position (1-based, inclusive) of the exon.  Must be > 0.
     * @param strand The {@link Strand} on which the variant and exon are read.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @param offsetIndelAdjustment An adjustment added to account for bases lost / gained in an Indel event.
     * @return A {@link String} representing the codon change for the splice site represented by the given parameters.
     */
    public static String createSpliceSiteCodonChange(final int variantStart,
                                                     final int exonNumber,
                                                     final int exonStart,
                                                     final int exonEnd,
                                                     final Strand strand,
                                                     final int offsetIndelAdjustment) {

        ParamUtils.isPositive(variantStart, "Genomic positions must be > 0.");
        ParamUtils.isPositive(exonNumber, "Exon number must be > 0.");
        ParamUtils.isPositive(exonStart, "Genomic positions must be > 0.");
        ParamUtils.isPositive(exonEnd, "Genomic positions must be > 0.");
        assertValidStrand( strand );

        char sign = '-';
        int offset = exonStart - variantStart;
        if ( Math.abs(offset) > Math.abs(variantStart - exonEnd)) {
            offset = variantStart - exonEnd;
            sign = '+';
        }
        offset = Math.abs(offset);

        if (strand == Strand.NEGATIVE) {
            if ( sign == '+' ) {
                sign = '-';
            }
            else {
                sign = '+';
            }
        }

        // Add our indel adjustment here:
        if ( sign == '+' ) {
            offset += offsetIndelAdjustment;
        }
        else {
            offset -= offsetIndelAdjustment;
        }

        // Make sure we correctly adjust for the zero crossing:
        if ( offset < 0 ) {
            offset *= -1;
            if ( sign == '+') {
                sign = '-';
            }
            else {
                sign = '+';
            }
        }

        return "c.e" + exonNumber + sign + offset;
    }

    /**
     * Creates the string representation of the codon change for the given {@link SequenceComparison}.
     * Requires that the given {@code seqComp} has the following fields defined with values that are not {@code null}:
     *     referenceAminoAcidSequence
     *     proteinChangeStartPosition
     *     proteinChangeEndPosition
     *     alternateAminoAcidSequence
     *     referenceAllele
     *     alternateAllele
     * In the case of an insertion, the given {@code seqComp} must have the following fields defined with values that are not {@code null}:
     *     codingSequenceAlleleStart
     *     transcriptCodingSequence
     *     strand  (must also not be {@link Strand#NONE}
     * @param seqComp {@link SequenceComparison} representing the alternate and reference alleles for a DNA sequence.  Must not be {@code null}.
     * @return A {@link String} representing the codon change for the given {@link SequenceComparison}.
     */
    public static String getProteinChangeString(final SequenceComparison seqComp) {

        Utils.nonNull(seqComp);
        Utils.nonNull(seqComp.getReferenceAminoAcidSequence());
        Utils.nonNull(seqComp.getProteinChangeStartPosition());
        Utils.nonNull(seqComp.getProteinChangeEndPosition());
        Utils.nonNull(seqComp.getAlternateAminoAcidSequence());
        Utils.nonNull(seqComp.getReferenceAllele());
        Utils.nonNull(seqComp.getAlternateAllele());

        if ( seqComp.getReferenceAllele().length() == 0 ) {
            throw new UserException("Reference allele cannot be empty.");
        }
        if ( seqComp.getAlternateAllele().length() == 0 ) {
            throw new UserException("Alternate allele cannot be empty.");
        }

        final String refAaSeq = seqComp.getReferenceAminoAcidSequence();
        final String altAaSeq = seqComp.getAlternateAminoAcidSequence();
        final Integer protChangeStartPos = seqComp.getProteinChangeStartPosition();
        final Integer protChangeEndPos = seqComp.getProteinChangeEndPosition();

        // Check for the ONP case:
        if ( GATKVariantContextUtils.isXnp(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {
            return getProteinChangeStringForOnp(refAaSeq, altAaSeq, protChangeStartPos, protChangeEndPos);
        }
        else {

            // Insertion:
            if ( GATKVariantContextUtils.isInsertion(seqComp.getAlignedReferenceAllele(), seqComp.getAlignedAlternateAllele()) ) {

                if ( GATKVariantContextUtils.isFrameshift(seqComp.getAlignedReferenceAllele(), seqComp.getAlignedAlternateAllele()) ) {

                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {
                        final String nextRefCodon = getNextReferenceCodon(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());
                        return "p." + createAminoAcidSequence(nextRefCodon) + (protChangeStartPos + 1) + "fs";
                    }
                    else {
                        // Right now, we don't care about anything that comes after the frameshift because it'll be garbage.
                        // Therefore, since we know that we're "inside" the effected codon, we report it as the frameshifted
                        // codon:
                        return "p." + createAminoAcidSequence(seqComp.getAlignedCodingSequenceReferenceAllele()) + protChangeStartPos + "fs";
                    }
                }
                else {
                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {
//                        final String nextRefCodon = getNextReferenceCodon(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());
//
//                         Here we adjust everything "right" by 3 bases (1 codon) because of the leading base that is
//                         required for indels:
//                        return "p." + protChangeStartPos + "_" + protChangeEndPos +
//                                createAminoAcidSequence(nextRefCodon) + ">" + createAminoAcidSequence(seqComp.getAlignedAlternateAllele().substring(3) + nextRefCodon);
                        return "p." + protChangeStartPos + "_" + (protChangeStartPos + 1) + "ins"+
                                createAminoAcidSequence(seqComp.getAlignedAlternateAllele().substring(3));
                    }
                    else {
                        // TODO: must do head / tail matching on the protein strings to deterimine what the positions of the protein change string actually are.
                        //       That is, you can have a mutation in protein position 979 that results in an amino acid insertion
                        //       between pp978 and pp979:
                        //          g.chr3:178948163_178948164insTGA ==>> c.(2935-2937)agg>aTGAgg ==>> R->MR ==>> p.978_979insM

                        if ( protChangeStartPos.equals(protChangeEndPos) ) {
                            // We have an insertion in a single base here.
                            // We have a different notation for this:
                            return "p." + protChangeStartPos + "ins" + createAminoAcidSequence(seqComp.getAlignedAlternateAllele());
                        }
                        else {

                            // Get the number of extra codons to get:
                            // NOTE: We only ever get 1 extra codon if we might overlap with the next codon
                            //       (as per oncotator conventions):
                            final int numAdditionalCodonsToGet = ((seqComp.getCodingSequenceAlleleStart() + 1) % 3 == 0 ? 1 : 0);

                            // Get the next few required reference codons:
                            final List<String> nextRefCodons = getNextReferenceCodons(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand(), numAdditionalCodonsToGet);


                            return "p." + protChangeStartPos + "_" + (protChangeStartPos + numAdditionalCodonsToGet) +
                                    createAminoAcidSequence(seqComp.getAlignedCodingSequenceReferenceAllele() + String.join("", nextRefCodons)) + ">" + createAminoAcidSequence(seqComp.getAlignedAlternateAllele() + String.join("", nextRefCodons));
                        }
                    }
                }

            }
            // Deletion:
            else {
                if ( GATKVariantContextUtils.isFrameshift(seqComp.getAlignedReferenceAllele(), seqComp.getAlignedAlternateAllele()) ) {
                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {

                        // Check to see if the deletion actually starts in the next codon, if it does then we skip the
                        // current codon and only display the next one:
                        if ((seqComp.getCodingSequenceAlleleStart() % 3) == 0) {
                            final String nextRefCodon = getNextReferenceCodon(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());
                            return "p." + createAminoAcidSequence(nextRefCodon) + (protChangeStartPos + 1) + "fs";
                        }
                        else {
                            return "p." + createAminoAcidSequence(seqComp.getAlignedCodingSequenceReferenceAllele()) + protChangeStartPos + "fs";
                        }
                    }
                    else {

                        // We started in the middle of a codon.
                        // We should report the current codon as the start position.

                        // Check to see if the deletion actually starts in the next codon, if it does then we skip the
                        // current codon and only display the next one:
                        if ((seqComp.getCodingSequenceAlleleStart() % 3) == 0) {
                            // First base (i.e. the required base before indels) is the last base of the first reference
                            // codon.  So we can skip it.
                            return "p." + refAaSeq.substring(1) + (protChangeStartPos + 1) + "fs";
                        }
                        else {
                            // Get the offset for where our REAL difference occurs:
                            final int refAaOffset = getOffsetForRefAaSequence(refAaSeq, altAaSeq);

                            // Return the string:
                            return "p." + refAaSeq.substring(refAaOffset) + (protChangeStartPos + refAaOffset) + "fs";
                        }
                    }
                }
                else {

                    // Determine how many codons to get:
                    final int numAdditionalCodonsToGet = (int) Math.ceil((seqComp.getAlternateAllele().length() + ((seqComp.getAlignedCodingSequenceAlleleStart() - seqComp.getCodingSequenceAlleleStart()) % 3)) / 3);

                    // Get the next few required reference codons:
                    final List<String> nextRefCodons = getNextReferenceCodons(seqComp.getTranscriptCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand(), numAdditionalCodonsToGet);

                    // Get the end position of the variant:
                    final int endPos = seqComp.getAlignedReferenceAlleleStop() + (3 * numAdditionalCodonsToGet);

                    // Get the string for all reference codons (including the next few codons we need):
                    final String allRefCodons = (seqComp.getAlignedCodingSequenceReferenceAllele() + String.join("", nextRefCodons)).toLowerCase();

                    // This means that the deletion is aligned with a codon:
                    if ( isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getReferenceAllele()) ) {

                        // An entire codon / set of codons got deleted.
                        // Just put `DEL` after it.
                        // NOTE: We need to ignore the first codon in the list because of the leading base for indels
//                        return "p." + (protChangeStartPos + 1) + "_" + endPos
//                                + createAminoAcidSequence(allRefCodons.substring(3)) + "del";
                        return "p." + createAminoAcidSequence(allRefCodons.substring(3)) + (protChangeStartPos + 1) + "del";
                    }
                    // Non-frameshift deletion starting within a codon:
                    else {
                        // We must report all old codons and then the new codons.
                        // However, for the new codons, we need to create a codon string that includes the new alternate
                        // allele and ONLY the bases that follow it (rounding to the next codon boundary).

                        final String referenceAminoAcidString = createAminoAcidSequence(allRefCodons);
                        final String alternateAminoAcidString = createAminoAcidSequence(seqComp.getAlignedCodingSequenceAlternateAllele() + String.join("", nextRefCodons));

                        final StringBuilder refSb = new StringBuilder();
                        final StringBuilder altSb = new StringBuilder();

                        final int referenceOffset = collapseNonFrameshiftProteinChangeStrings( referenceAminoAcidString,
                                                                                alternateAminoAcidString,
                                                                                refSb,
                                                                                altSb );

                        // Use the offset to adjust the starting positions:
                        int adjustedProteinStartPos = protChangeStartPos;
                        int adjustedProteinEndPos = protChangeEndPos;
                        if ( referenceOffset > 0) {
                            adjustedProteinStartPos += referenceOffset;
                        }
                        if ( referenceOffset < 0) {
                            adjustedProteinEndPos += referenceOffset;
                        }

                        // Check for "pure" deletion:
                        if ( altSb.length() == 0 ) {
                            return "p." + refSb.toString() + adjustedProteinStartPos + "del";
                        }
                        else {
                            return "p." + adjustedProteinStartPos + "_" + adjustedProteinEndPos
                                    + refSb.toString() + ">" + altSb.toString();
                        }
                    }
                }
            }
        }
    }

    private static int collapseNonFrameshiftProteinChangeStrings( final String refAaSeq,
                                                                  final String altAaSeq,
                                                                  final StringBuilder outRefAa,
                                                                  final StringBuilder outAltAa) {
        ParamUtils.isPositive( refAaSeq.length(), "Reference amino acid string must be of non-zero length." );
        ParamUtils.isPositive( altAaSeq.length(), "Alternate amino acid string must be of non-zero length." );

        int offset = 0;

        // Check for front overlap:
        if ( refAaSeq.charAt(0) == altAaSeq.charAt(0) ) {

            while ( (offset < refAaSeq.length()) && (offset < altAaSeq.length()) &&
                    (refAaSeq.charAt(offset) == altAaSeq.charAt(offset)) ) {
                ++offset;
            }

            for ( int i = offset; i < refAaSeq.length() ; ++i ) {
                outRefAa.append( refAaSeq.charAt(i) );
            }

            for ( int i = offset; i < altAaSeq.length() ; ++i ) {
                outAltAa.append( altAaSeq.charAt(i) );
            }

        }
        // Check for back overlap:
        else if ( refAaSeq.charAt(refAaSeq.length()-1) == altAaSeq.charAt(altAaSeq.length()-1) ) {

            while ( (offset < refAaSeq.length()) && (offset < altAaSeq.length()) &&
                    (refAaSeq.charAt(refAaSeq.length() - 1 - offset) == altAaSeq.charAt(altAaSeq.length()  - 1 - offset)) ) {
                ++offset;
            }

            for ( int i = 0; i < (refAaSeq.length() - offset) ; ++i ) {
                outRefAa.append( refAaSeq.charAt(i) );
            }

            for ( int i = 0; i < (altAaSeq.length() - offset) ; ++i ) {
                outAltAa.append( altAaSeq.charAt(i) );
            }
            offset = -offset;
        }

        return offset;
    }

    private static String collapseNonFrameshiftProteinChangeString(final int proteinChangeStartPosition, final int proteinChangeEndPosition,
                                                                   final String refAaSeq, final String altAaSeq ) {
        if ( proteinChangeStartPosition == proteinChangeEndPosition ) {
            if ( altAaSeq.isEmpty() ) {
                return "p." + refAaSeq + proteinChangeStartPosition + "del";
            }
            else {
                return "p." + proteinChangeStartPosition + "_" + proteinChangeEndPosition + "ins"
                        + altAaSeq;
            }
        }
        else {
            return "p." + proteinChangeStartPosition + "_" + proteinChangeEndPosition
                    + refAaSeq + ">" + altAaSeq;
        }
    }

    /**
     * Get the offset of the first different base in the given reference amino acid sequence.
     * @param refAaSeq {@link String} containing the single-letter abbreviations for the amino acids in the reference allele's coding sequence.
     * @param altAaSeq {@link String} containing the single-letter abbreviations for the amino acids in the alternate allele's coding sequence.
     * @return The position of the first base in {@code refAaSeq} that is not the same in altAaSeq.
     */
    private static int getOffsetForRefAaSequence(final String refAaSeq, final String altAaSeq) {

        Utils.nonNull(refAaSeq);
        Utils.nonNull(altAaSeq);

        ParamUtils.isPositive(refAaSeq.length(), "Reference amino acid sequence must be of length greater than 0!");

        // Degenerate cases:
        if ( altAaSeq.isEmpty() ) {
            return 0;
        }
        else {
            int index = 0;
            while ( (index < refAaSeq.length()) && (index < altAaSeq.length()) ) {
                if ( refAaSeq.charAt(index) != altAaSeq.charAt(index) ) {
                    break;
                }
                ++index;
            }
            return index;
        }
    }

    /**
     * Gets the {@link String} representing the protein change for an ONP.
     * @param referenceAminoAcidSequence A {@link String} containing the single-letter codes for the {@link AminoAcid}s in the Reference Allele.
     * @param alternateAminoAcidSequence A {@link String} containing the single-letter codes for the {@link AminoAcid}s in the Alternate Allele.
     * @param protChangeStartPos The position (1-based, inclusive) of the start of the protein change in the sequence of proteins coded by the transcript containing the given alleles.
     * @param protChangeEndPos The position (1-based, inclusive) of the start of the protein change in the sequence of proteins coded by the transcript containing the given alleles.
     * @return The {@link String} representing the protein change for the given alleles / amino acid sequences.
     */
    private static String getProteinChangeStringForOnp(final String referenceAminoAcidSequence,
                                                       final String alternateAminoAcidSequence,
                                                       final int protChangeStartPos,
                                                       final int protChangeEndPos) {

        // We should go through our strings and make sure we only render the parts of the protein that have actually
        // changed.  This means that we keep track of the `same` and `different` parts of the sequences.
        // `Same` parts at the front & back of the string are ignored.
        // We keep them in the middle.
        // Then we recompute the position in which the protein has changed.

        boolean foundStartDiff = false;
        boolean foundEndDiff = false;

        int startingDifference = 0;
        int endingDifference = 0;

        int startPos = protChangeStartPos;
        int endPos   = protChangeEndPos;

        for ( int i = 0 ; i < referenceAminoAcidSequence.length() ; ++i ) {
            final int rIndx = referenceAminoAcidSequence.length() - 1 - i;
            if ( (!foundStartDiff) && (referenceAminoAcidSequence.charAt(i) != alternateAminoAcidSequence.charAt(i)) ) {
                startingDifference = i;
                foundStartDiff = true;
            }
            if ( (!foundEndDiff) && (referenceAminoAcidSequence.charAt(rIndx) != alternateAminoAcidSequence.charAt(rIndx)) ) {
                endingDifference = i;
                foundEndDiff = true;
            }
        }

        // Set the new start / stop positions:
        startPos += startingDifference;
        endPos -= endingDifference;

        // Set the new ref and alt amino acid sequences:
        final String refAaSeq = referenceAminoAcidSequence.substring(startingDifference, referenceAminoAcidSequence.length() - endingDifference);
        final String altAaSeq = alternateAminoAcidSequence.substring(startingDifference, alternateAminoAcidSequence.length() - endingDifference);

        if (startPos == endPos) {
            return "p." + refAaSeq + startPos +
                    altAaSeq;
        } else {
            return "p." + startPos
                    + "_" + endPos + refAaSeq + '>' + altAaSeq;
        }
    }

    /**
     * Gets the {@link String} representing the protein change for a Deletion.
     * @param referenceAllele The {@link String} of bases contained in the Reference allele.
     * @param alternateAllele The {@link String} of bases contained in the Alternate allele.
     * @param referenceAminoAcidSequence A {@link String} containing the single-letter codes for the {@link AminoAcid}s in the Reference Allele.
     * @param protChangeStartPos The position (1-based, inclusive) of the start of the protein change in the sequence of proteins coded by the transcript containing the given alleles.
     * @return The {@link String} representing the protein change for the given alleles / amino acid sequences.
     */
    private static String getProteinChangeStringForDeletion(final String referenceAllele,
                                                            final String alternateAllele,
                                                            final String referenceAminoAcidSequence,
                                                            final int protChangeStartPos) {
        // NOTE: We know that because of the way deletions are created
        // (i.e. a base is added to the beginning of the string), we can ignore the first amino acid (3 bases) of the
        // reference:
        String decorator = "del";
        if ( GATKVariantContextUtils.isFrameshift(referenceAllele, alternateAllele) ) {
            decorator = "fs";
        }
        return "p." + referenceAminoAcidSequence.substring(1) + (protChangeStartPos + 1) + decorator;
    }


    /**
     * Gets the {@link String} representing the protein change for a Frameshift Insertion.
     * @param codingSequenceAlleleStart The position (1-based, inclusive) of the start of the coding sequence of the alleles.  Not guaranteed to be in-frame.
     * @param refAaSeq A {@link String} containing the single-letter codes for the {@link AminoAcid}s in the Reference Allele.
     * @param protChangeStartPos The position (1-based, inclusive) of the start of the protein change in the sequence of proteins coded by the transcript containing the given alleles.
     * @param nextRefAaSeq The {@link String} containing the next single-letter codes for the {@link AminoAcid}s in the Reference Allele
     * @return The {@link String} representing the protein change for the given alleles / amino acid sequences.
     */
    private static String getProteinChangeStringForInsertionFrameshift(final int codingSequenceAlleleStart,
                                                                       final String refAaSeq,
                                                                       final int protChangeStartPos,
                                                                       final String nextRefAaSeq) {
        // If the variant position is in frame, then we have to use the next amino acid, because it will be the first
        // amino acid that is affected.
        // We add 1 because of the convention that the variant occurs just after the base specified for insertions.
        // TODO: This is really bad - we are tying our output to a specific input format.  FIX IT.
        if ( isPositionInFrame(codingSequenceAlleleStart + 1) ) {
            return "p." + nextRefAaSeq + (protChangeStartPos + 1) + "fs";
        }
        else {
            // Out of frame insertions affect the current protein position, so the current position should be
            // used.
            return "p." + refAaSeq + protChangeStartPos + "fs";
        }
    }

    /**
     * Gets the {@link String} representing the protein change for an in-frame Insertion with an in-frame start position.
     * @param protChangeStartPos The position (1-based, inclusive) of the start of the protein change in the sequence of proteins coded by the transcript containing the given alleles.
     * @param altAaSeq A {@link String} containing the single-letter codes for the {@link AminoAcid}s in the Alternate Allele.
     * @return The {@link String} representing the protein change for the given alleles / amino acid sequences.
     */
    private static String getProteinChangeStringForInsertionInFrameWithInFrameStartPosition(final int protChangeStartPos,
                                                                                            final String altAaSeq) {
        // NOTE: We know that because of the way insertions are created
        // (i.e. a base is added to the beginning of the string), we can ignore the first amino acid (3 bases) of the
        // amino acid sequence:
        return "p." + protChangeStartPos + "_" + (protChangeStartPos+1) + "ins" + altAaSeq.substring(1);
    }

    /**
     * Gets the {@link String} representing the protein change for an in-frame Insertion with an out-of-frame start position.
     * @param refAaSeq A {@link String} containing the single-letter codes for the {@link AminoAcid}s in the Reference Allele.
     * @param altAaSeq A {@link String} containing the single-letter codes for the {@link AminoAcid}s in the Alternate Allele.
     * @param protChangeStartPos The position (1-based, inclusive) of the start of the protein change in the sequence of proteins coded by the transcript containing the given alleles.
     * @param alternateAlleleLength The length (in bases) of the raw Alternate allele as given by the input file.
     * @return The {@link String} representing the protein change for the given alleles / amino acid sequences.
     */
    private static String getProteinChangeStringForInsertionInFrameWithOutOfFrameStartPosition(final String refAaSeq,
                                                                                               final String altAaSeq,
                                                                                               final int protChangeStartPos,
                                                                                               final int alternateAlleleLength ) {
        // For this case we know that protChangeStartPos is the position of the current amino acid in which
        // the insertion took place.
        // We also know that refAaSeq is the amino acid at protChangeStartPos.


        // To get the protein change itself, we need to iterate through the ref and alt amino acid sequences
        // and determine which are overlaps (starting from the back).  Any overlaps are discarded until
        // differences are found.
        // Then the rest of the sequence is said to be new.
        final StringBuilder altAminoAcidBuilder = new StringBuilder();
        int i = refAaSeq.length()-1;
        int j = altAaSeq.length()-1;
        int skipped = 0;
        final int expectedLength = ((alternateAlleleLength-1) / 3);
        while ( altAminoAcidBuilder.length() < expectedLength ) {

            if ( (i < 0) || (refAaSeq.charAt(i) != altAaSeq.charAt(j)) ) {
                // Get all remaining amino acids into altAminoAcidBuilder:
                while ((j >= 0) && (altAminoAcidBuilder.length() != expectedLength)) {
                    altAminoAcidBuilder.append(altAaSeq.charAt(j--));
                }
            } else {
                --i;
                --j;
                ++skipped;
            }
        }
        altAminoAcidBuilder.reverse();

        // NOTE: We know that because of the way insertions are created
        // (i.e. a base is added to the beginning of the string), we can ignore the first amino acid (3 bases) of the
        // amino acid sequence:
        return "p." + (protChangeStartPos - skipped) + "_" + (protChangeStartPos - skipped + 1) + "ins" + altAminoAcidBuilder.toString();

    }

    /**
     * Get the coding sequence change string from the given {@link SequenceComparison}
     * This method is assumed to be called only when the variant occurs in a coding region of the genome.
     * Requires that the given {@code seqComp} has the following fields defined with values that are not {@code null}:
     *     codingSequenceAlleleStart
     *     referenceAllele
     *     alternateAllele
     * @param seqComp {@link SequenceComparison} from which to construct the coding sequence change string.  Must not be {@code null}.
     * @return A {@link String} representing the coding sequence change between the ref and alt alleles in {@code seqComp}.
     */
    public static String getCodingSequenceChangeString( final SequenceComparison seqComp ) {

        Utils.nonNull(seqComp);
        Utils.nonNull(seqComp.getCodingSequenceAlleleStart());
        Utils.nonNull(seqComp.getReferenceAllele());
        Utils.nonNull(seqComp.getAlternateAllele());

        // Check for ONP:
        if ( GATKVariantContextUtils.isXnp(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {
            if (seqComp.getAlternateAllele().length() > 1) {
                return "c." + seqComp.getCodingSequenceAlleleStart() + "_" + (seqComp.getCodingSequenceAlleleStart() + seqComp.getReferenceAllele().length() - 1) +
                        seqComp.getReferenceAllele() + ">" + seqComp.getAlternateAllele();
            } else {
                return "c." + seqComp.getCodingSequenceAlleleStart() +
                        seqComp.getReferenceAllele() + ">" + seqComp.getAlternateAllele();
            }
        }
        // Check for Insertion:
        else if ( GATKVariantContextUtils.isInsertion(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {
            return "c." + seqComp.getCodingSequenceAlleleStart() + "_" + (seqComp.getCodingSequenceAlleleStart() + 1) +
                    "ins" + seqComp.getAlternateAllele().substring(seqComp.getReferenceAllele().length()).toUpperCase();
        }
        // Must be a Deletion:
        else {

            int start = seqComp.getCodingSequenceAlleleStart() + seqComp.getAlternateAllele().length();
            int end = seqComp.getCodingSequenceAlleleStart() + seqComp.getReferenceAllele().length() - 1;

            // If we have exon information, and we SHOULD, we use it to trim the start/stop coordinates
            // of the cDNA string to the extants of the coding region:
            if ( (seqComp.getExonStartPosition() != null) && (seqComp.getExonEndPosition() != null) ) {
                final int cdsExonStart = seqComp.getCodingSequenceAlleleStart() - (seqComp.getAlleleStart() - seqComp.getExonStartPosition());
                final int cdsExonEnd   = cdsExonStart + (seqComp.getExonEndPosition() - seqComp.getExonStartPosition());

                if ( start < cdsExonStart ) {
                    start = cdsExonStart;
                }
                if ( end > cdsExonEnd ) {
                    end = cdsExonEnd;
                }
            }

            if ( start == end ) {
                return "c." + start + "del" +
                        seqComp.getReferenceAllele().substring(seqComp.getAlternateAllele().length()).toUpperCase();
            }
            else {
                return "c." + start + "_" + end + "del" +
                        seqComp.getReferenceAllele().substring(seqComp.getAlternateAllele().length()).toUpperCase();
            }
        }
    }

    /**
     * Get the coding sequence change string from the given {@link SequenceComparison}
     * @param transcriptPosition The position in the transcript of the given splice site variant.  Must be > 0.
     * @return A {@link String} representing the coding sequence change between the ref and alt alleles in {@code seqComp}.
     */
    public static String getCodingSequenceChangeStringForExonSpliceSite( final int transcriptPosition ) {

        ParamUtils.isPositive( transcriptPosition, "Transcript position must be > 0." );

        if ( transcriptPosition < 1 ) {
            throw new GATKException("Encountered transcript position less than 1 (transcript positions are 1-based): " + transcriptPosition + " < " + 1);
        }

        return "c." + transcriptPosition + "_splice";
    }

    /**
     * Creates an amino acid sequence from a given coding sequence.
     * If the coding sequence is not evenly divisible by 3, the remainder bases will not be included in the coding sequence.
     * @param codingSequence The coding sequence from which to create an amino acid sequence.  Must not be {@code null}.
     * @return A {@link String} containing a sequence of single-letter amino acids.
     */
    public static String createAminoAcidSequence(final String codingSequence) {

        Utils.nonNull(codingSequence);

        final StringBuilder sb = new StringBuilder();

        // Ensure that we don't have remainder bases:
        int maxIndex = codingSequence.length();
        if ( maxIndex % 3 != 0 ) {
            maxIndex = (int)Math.floor(maxIndex / 3) * 3;
            logger.warn("createAminoAcidSequence given a coding sequence of length not divisible by 3.  Dropping bases from the end: " + (codingSequence.length() % 3));
        }

        for ( int i = 0; i < maxIndex; i += 3 ) {
            final AminoAcid aa = getEukaryoticAminoAcidByCodon(codingSequence.substring(i, i+3));
            if ( aa == null ) {
                throw new UserException.MalformedFile("File contains a bad codon sequence that has no amino acid equivalent: " + codingSequence.substring(i, i+3));
            }
            else {
                sb.append(aa.getLetter());
            }
        }
        return sb.toString();
    }

    /**
     * Get the coding sequence-aligned allele based on stop and start position.
     * @param codingSequence Coding sequence from which the allele should be derived.  Must not be {@code null}.
     * @param alignedAlleleStart Start position of the allele (1-indexed, inclusive).  Must not be {@code null}.  Must be > 0.
     * @param alignedAlleleStop Stop position of the allele (1-indexed, inclusive).  Must not be {@code null}.  Must be > 0.
     * @param strand {@link Strand} on which the allele is coded.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return The {@link String} representation of the allele.
     */
    private static String getAlignedAlleleSequence(final String codingSequence,
                                                   final Integer alignedAlleleStart,
                                                   final Integer alignedAlleleStop,
                                                   final Strand strand) {
        Utils.nonNull(codingSequence);
        Utils.nonNull(alignedAlleleStart);
        ParamUtils.isPositive( alignedAlleleStart, "Genome positions must be > 0." );
        Utils.nonNull(alignedAlleleStop);
        ParamUtils.isPositive( alignedAlleleStop, "Genome positions must be > 0." );
        assertValidStrand( strand );

        // Get our indices:
        // Subtract 1 because we're 1-based.
        int start = alignedAlleleStart - 1;
        int end = alignedAlleleStop;

        final String alignedAlleleSeq;

        if ( strand == Strand.POSITIVE ) {

            if ( end > codingSequence.length() ) {
                throw new TranscriptCodingSequenceException("Gencode transcript ends at position " + end + " but codingSequence is only " + codingSequence.length() + " bases long!");
            }
            else {
                alignedAlleleSeq = codingSequence.substring(start, end);
            }
        }
        else {
            // Negative strand means we need to reverse complement and go from the other end:
            start = codingSequence.length() - alignedAlleleStop;
            end = codingSequence.length() - alignedAlleleStart + 1;

            if ( end > codingSequence.length() ) {
                throw new TranscriptCodingSequenceException("Gencode transcript ends at position " + end + " but codingSequence is only " + codingSequence.length() + " bases long!");
            }
            else {
                alignedAlleleSeq = ReadUtils.getBasesReverseComplement(codingSequence.substring(start, end).getBytes());
            }
        }

        return alignedAlleleSeq;
    }

    /**
     * Gets the coding sequence for the allele with given start and stop positions, codon-aligned to the start of the reference sequence.
     * @param codingSequence The whole coding sequence for this transcript.  Must not be {@code null}.
     * @param alignedAlleleStart The codon-aligned position (1-based, inclusive) of the allele start.  Must not be {@code null}.
     * @param alignedAlleleStop The codon-aligned position (1-based, inclusive) of the allele stop.  Must not be {@code null}.
     * @param refAllele The reference {@link Allele}.  Must not be {@code null}.
     * @param refAlleleStart The position (1-based, inclusive) of where the reference allele starts.  Must not be {@code null}.  Must be > 0.
     * @param strand The {@link Strand} on which the alleles are found.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return A {@link String} containing the reference allele coding sequence.
     */
    public static String getAlignedCodingSequenceAllele(final String codingSequence,
                                                        final Integer alignedAlleleStart,
                                                        final Integer alignedAlleleStop,
                                                        final Allele  refAllele,
                                                        final Integer refAlleleStart,
                                                        final Strand strand) {
        Utils.nonNull(codingSequence);
        Utils.nonNull(alignedAlleleStart);
        ParamUtils.isPositive( alignedAlleleStart, "Genome positions must be > 0." );
        Utils.nonNull(alignedAlleleStop);
        ParamUtils.isPositive( alignedAlleleStop, "Genome positions must be > 0." );
        Utils.nonNull(refAllele);
        Utils.nonNull(refAlleleStart);
        ParamUtils.isPositive( refAlleleStart, "Genome positions must be > 0." );

        assertValidStrand( strand );

        String alignedAlleleSeq = getAlignedAlleleSequence(codingSequence, alignedAlleleStart, alignedAlleleStop, strand);

        // Check whether our reference sequence is derived from the reference or if it should be derived from the given
        // reference.
        final String expectedReferenceSequence;
        if ( strand == Strand.POSITIVE ) {
            expectedReferenceSequence = codingSequence.substring(refAlleleStart - 1, refAlleleStart - 1 + refAllele.length());
        }
        else {
            final int start = codingSequence.length() - (refAlleleStart - 1 + refAllele.length());
            final int end = codingSequence.length() - refAlleleStart;
            expectedReferenceSequence = ReadUtils.getBasesReverseComplement( codingSequence.substring(start, end).getBytes() );
        }

        // NOTE: This check appears to be reduntant, but in actuality, it is required.
        //       Because we reconstruct the coding sequence allele separately from the reference allele, we need to check this
        //       again to make sure we have the right alleles given our input.
        if ( !expectedReferenceSequence.equals(refAllele.getBaseString()) ) {
            // Oh noes!
            // Ref allele is different from reference sequence!
            // Oh well, we should use the reference we were given anyways...
            final String substitutedAlignedSeq = getAlternateSequence(codingSequence, refAlleleStart, refAllele, refAllele);

            // We use the positive strand here because we have already reverse complemented the sequence in the call
            // above.
            final String substitutedAlignedAlleleSeq = getAlignedAlleleSequence(substitutedAlignedSeq, alignedAlleleStart, alignedAlleleStop, Strand.POSITIVE);

            // Warn the user!
            logger.warn("Reference allele is different than the reference coding sequence!  Substituting given allele for sequence code (" + alignedAlleleSeq + "->" + substitutedAlignedAlleleSeq + ")");

            // Set up our return value:
            alignedAlleleSeq = substitutedAlignedAlleleSeq;
        }

        return alignedAlleleSeq;
    }

    /**
     * Get the aligned coding sequence for the given reference allele.
     * @param referenceSnippet {@link String} containing a short excerpt of the reference sequence around the given reference allele.  Must not be {@code null}.
     * @param referencePadding Number of bases in {@code referenceSnippet} before the reference allele starts.  This padding exists at the end as well (plus some other bases to account for the length of the alternate allele if it is longer than the reference).  Must be >= 0.
     * @param refAllele The reference {@link Allele}.  Must not be {@code null}.
     * @param codingSequenceRefAlleleStart The position (1-based, inclusive) in the coding sequence where the {@code refAllele} starts.  Must be > 0.
     * @param alignedRefAlleleStart The in-frame position (1-based, inclusive) of the first base of the codon containing the reference allele.  Must be > 0.
     * @return A {@link String} of in-frame codons that contain the entire reference allele.
     */
    public static String getAlignedRefAllele(final String referenceSnippet,
                                             final int referencePadding,
                                             final Allele refAllele,
                                             final int codingSequenceRefAlleleStart,
                                             final int alignedRefAlleleStart ) {

        Utils.nonNull(referenceSnippet);
        Utils.nonNull(refAllele);

        ParamUtils.isPositiveOrZero( referencePadding, "Padding must be >= 0." );
        ParamUtils.isPositive( codingSequenceRefAlleleStart, "Genome positions must be > 0." );
        ParamUtils.isPositive( alignedRefAlleleStart, "Genome positions must be > 0." );

        final int extraBasesNeeded = (codingSequenceRefAlleleStart - alignedRefAlleleStart);
        int refStartPos = referencePadding - extraBasesNeeded;

        if ( refStartPos < 0 ) {
            refStartPos = 0;
        }

        // Round to the nearest multiple of 3 to get the end position.
        int refEndPos = refStartPos + (int)(Math.ceil((extraBasesNeeded + refAllele.length()) / 3.0) * 3);

        // Create the aligned reference:
        String alignedReferenceAllele = referenceSnippet.substring(refStartPos, refEndPos);

        // Make sure our reference is what we expect it to be with the ref allele:
        if ( !alignedReferenceAllele.substring(extraBasesNeeded, extraBasesNeeded + refAllele.length()).equals(refAllele.getBaseString()) ) {
            // Oh noes!
            // Ref allele is different from reference sequence!

            // Oh well, we should use the reference we were given anyways...
            final String substitutedReferenceSnippet = getAlternateSequence(referenceSnippet, referencePadding + 1, refAllele, refAllele);
            refEndPos = refStartPos + (int)(Math.ceil((extraBasesNeeded + refAllele.length()) / 3.0) * 3);

            final String substitutedAlignedAlleleSeq = substitutedReferenceSnippet.substring(refStartPos, refEndPos);

            // Warn the user!
            logger.warn("Reference allele is different than the reference coding sequence!  Substituting given allele for sequence code (" + alignedReferenceAllele + "->" + substitutedAlignedAlleleSeq + ")");

            // Set up our return value:
            alignedReferenceAllele = substitutedAlignedAlleleSeq;
        }

        return alignedReferenceAllele;

    }

    /**
     * Get a string of bases around a variant (specified by reference and alternate alleles), including the reference allele itself.
     * ASSUMES: that the given {@link ReferenceContext} is already centered on the variant location.
     * @param refAllele The reference {@link Allele} for the variant in question.  If on {@link Strand#NEGATIVE}, must have already been reverse complemented.  Must not be {@code null}.
     * @param altAllele The alternate {@link Allele} for the variant in question.  If on {@link Strand#NEGATIVE}, must have already been reverse complemented.  Must not be {@code null}.
     * @param strand The {@link Strand} on which the variant in question lives.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @param referenceWindowInBases The number of bases to the left and right of the variant to return.  Must be > 0.
     * @param referenceContext The {@link ReferenceContext} centered around the variant in question.  Must not be {@code null}.
     * @return A string containing {@code referenceWindowInBases} bases to either side of the specified refAllele.
     */
    public static String getBasesInWindowAroundReferenceAllele( final Allele refAllele,
                                                                final Allele altAllele,
                                                                final Strand strand,
                                                                final int referenceWindowInBases,
                                                                final ReferenceContext referenceContext) {
        Utils.nonNull( refAllele );
        Utils.nonNull( altAllele );
        assertValidStrand( strand );
        Utils.nonNull( referenceContext );

        // Calculate our window to include any extra bases but also have the right referenceWindowInBases:
        final int endWindow = refAllele.length() >= altAllele.length() ? referenceWindowInBases + refAllele.length() - 1 : referenceWindowInBases + altAllele.length() - 1;

        final String referenceBases;

        if ( strand == Strand.POSITIVE ) {
            // Get the reference sequence:
            referenceBases = new String(referenceContext.getBases(referenceWindowInBases, endWindow));
        }
        else {
            // Get the reference sequence:
            referenceBases = ReadUtils.getBasesReverseComplement(referenceContext.getBases(referenceWindowInBases, endWindow));
        }

        return referenceBases;
    }

    /**
     * Get the Protein change start position (1-based, inclusive) given the aligned position of the coding sequence.
     * @param alignedCodingSequenceAlleleStart Position (1-based, inclusive) of the start of the allele in the coding sequence.  Must not be {@code null}.  Must be > 0.
     * @return The position (1-based, inclusive) of the protein change in the amino acid sequence.
     */
    public static int getProteinChangePosition(final Integer alignedCodingSequenceAlleleStart) {

        Utils.nonNull(alignedCodingSequenceAlleleStart);
        ParamUtils.isPositive( alignedCodingSequenceAlleleStart, "Genome positions must be > 0." );
        return ((alignedCodingSequenceAlleleStart-1) / 3) + 1; // Add 1 because we're 1-based.
    }

    /**
     * Get the Protein change end position (1-based, inclusive) given the protein change start position and aligned alternate allele length.
     * @param proteinChangeStartPosition Position (1-based, inclusive) of the start of the protein change in the amino acid sequence.  Must not be {@code null}.  Must be > 0.
     * @param alignedAlternateAlleleLength Length of the aligned alternate allele in bases.  Must not be {@code null}.  Must be > 0.
     * @return The position (1-based, inclusive) of the end of the protein change in the amino acid sequence.
     */
    public static int getProteinChangeEndPosition(final Integer proteinChangeStartPosition, final Integer alignedAlternateAlleleLength) {

        Utils.nonNull(proteinChangeStartPosition);
        Utils.nonNull(alignedAlternateAlleleLength);

        ParamUtils.isPositive( proteinChangeStartPosition, "Genome positions must be > 0." );
        ParamUtils.isPositive( alignedAlternateAlleleLength, "Allele length > 0." );

        // We subtract 1 because we're 1-based.
        return proteinChangeStartPosition + getProteinChangePosition(alignedAlternateAlleleLength) - 1;
    }

    /**
     * Get the full alternate sequence given a reference coding sequence, and two alleles.
     * @param referenceSequence The reference sequence on which to base the resulting alternate sequence.  Must not be {@code null}.
     * @param alleleStartPos Starting position (1-based, inclusive) for the ref and alt alleles in the given {@code referenceSequence}.  Must be > 0.
     * @param refAllele Reference Allele.  Used for the length of the reference (content ignored).  Must not be {@code null}.
     * @param altAllele Alternate Allele.  Used for both content and length of the alternate allele.  Must not be {@code null}.
     * @return The coding sequence that includes the given alternate allele in place of the given reference allele.
     */
    public static String getAlternateSequence(final String referenceSequence,
                                              final int alleleStartPos,
                                              final Allele refAllele,
                                              final Allele altAllele ) {

        Utils.nonNull(referenceSequence);
        Utils.nonNull(refAllele);
        Utils.nonNull(altAllele);

        ParamUtils.isPositive( alleleStartPos, "Genome positions must be > 0." );

        // We have to subtract 1 here because we need to account for the 1-based indexing of
        // the start and end of the coding region:
        final int alleleIndex = Math.abs(alleleStartPos - 1);

        return referenceSequence.substring(0, alleleIndex) +
                altAllele.getBaseString() +
                referenceSequence.substring(alleleIndex + refAllele.length());
    }

    /**
     * Gets the start position relative to the start of the coding sequence for a variant based on the given {@code variantGenomicStartPosition}.
     * It is assumed:
     *      {@code codingRegionGenomicStartPosition} <= {@code variantGenomicStartPosition} <= {@code codingRegionGenomicEndPosition}
     * The transcript start position is the genomic position the transcript starts assuming `+` traversal.  That is, it is the lesser of start position and end position.
     * This is important because we determine the relative positions based on which direction the transcript is read.
     * @param variantGenomicStartPosition Start position (1-based, inclusive) of a variant in the genome.  Must be > 0.
     * @param codingRegionGenomicStartPosition Start position (1-based, inclusive) of a transcript in the genome.  Must be > 0.
     * @param codingRegionGenomicEndPosition End position (1-based, inclusive) of a transcript in the genome.  Must be > 0.
     * @param strand {@link Strand} from which strand the associated transcript is read.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return The start position (1-based, inclusive) relative to the start of the coding sequence of a variant.
     */
    public static int getTranscriptAlleleStartPosition(final int variantGenomicStartPosition,
                                                       final int codingRegionGenomicStartPosition,
                                                       final int codingRegionGenomicEndPosition,
                                                       final Strand strand) {

        ParamUtils.isPositive( variantGenomicStartPosition, "Genome positions must be > 0." );
        ParamUtils.isPositive( codingRegionGenomicStartPosition, "Genome positions must be > 0." );
        ParamUtils.isPositive( codingRegionGenomicEndPosition, "Genome positions must be > 0." );
        assertValidStrand( strand );

        if ( strand == Strand.POSITIVE ) {
            return variantGenomicStartPosition - codingRegionGenomicStartPosition + 1;
        }
        else {
            return codingRegionGenomicEndPosition - variantGenomicStartPosition + 1;
        }
    }

    /**
     * Creates and returns the coding sequence given a {@link ReferenceContext} and a {@link List} of {@link Locatable} representing a set of Exons.
     * Locatables start and end values are inclusive.
     * Assumes {@code exonList} ranges are indexed by 1.
     * @param reference A {@link ReferenceContext} from which to construct the coding region.  Must not be {@code null}.
     * @param exonList A {@link List} of {@link Locatable} representing a set of Exons to be concatenated together to create the coding sequence.  Each exon must have the same {@code contig} as {@code reference}.  Must not be {@code null}.
     * @param strand The {@link Strand} from which the exons are to be read.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return A string of bases for the given {@code exonList} concatenated together.
     */
    @Deprecated // Deprecated until fixed!
    public static String getCodingSequence(final ReferenceContext reference, final List<? extends Locatable> exonList, final Strand strand) {

        // TODO: Fix this - it's broken (and re-enable the tests).

        Utils.nonNull(reference);
        Utils.nonNull(exonList);
        assertValidStrand( strand );

        // Sanity check:
        if (exonList.size() == 0) {
            return "";
        }

        final StringBuilder sb = new StringBuilder();

        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;

        // Start by sorting our list of exons.
        // This is very important to ensure that we have all sequences in the right order at the end
        // and so we can support different read directions:
        exonList.sort((lhs, rhs) -> lhs.getStart() < rhs.getStart() ? -1 : (lhs.getStart() > rhs.getStart() ) ? 1 : 0 );

        for ( final Locatable exon : exonList ) {

            // First a basic sanity check:
            if ( !exon.getContig().equals(reference.getWindow().getContig()) ) {
                throw new GATKException("Cannot create a coding sequence! Contigs not the same - Ref: "
                        + reference.getInterval().getContig() + ", Exon: " + exon.getContig());
            }

            if ( start > exon.getStart() ) { start = exon.getStart(); }
            if ( end < exon.getEnd() ) { end = exon.getEnd(); }
        }

        // Get the window so we can convert to reference coordinates from genomic coordinates of the exons:
        final SimpleInterval refWindow = reference.getWindow();
        final byte[] bases = reference.getBases(Math.abs(start - reference.getInterval().getStart()), Math.abs(reference.getInterval().getEnd() - end));

        // If we're going in the opposite direction, we must reverse this reference base array
        // and complement it.
        if ( strand == Strand.NEGATIVE ) {
            for( int i = 0; i < bases.length / 2; i++) {
                final byte b = SequenceUtil.complement( bases[i] );
                bases[i] = SequenceUtil.complement( bases[bases.length - 1 - i] );
                bases[bases.length - 1 - i] = b;
            }
        }

        // Go through and grab our sequences based on our exons:
        for ( final Locatable exon : exonList ) {

            // Subtract 1 from start because positions are indexed by 1.
            int exonStartArrayCoord = exon.getStart() - refWindow.getStart() - 1;

            // Sanity check just in case the exon and ref window start at the same place:
            if ( exonStartArrayCoord == -1 ) {
                exonStartArrayCoord = 0;
            }

            // Add 1 to end because end range in copyOfRange is exclusive
            final int exonEndArrayCoord = exonStartArrayCoord + (exon.getEnd() - exon.getStart()) + 1;

            // TODO: find a better / faster way to do this:
            sb.append(
                    new String(
                            Arrays.copyOfRange(bases, exonStartArrayCoord, exonEndArrayCoord)
                    )
            );
        }

        return sb.toString();
    }

    /**
     * Get the overlapping exon start/stop as a {@link SimpleInterval} for the given altAllele / reference.
     * @param refAllele {@link Allele} for the given {@code altAllele}.  Must not be {@code null}.
     * @param altAllele {@link Allele} to locate on an exon.  Must not be {@code null}.
     * @param contig Contig on which the altAllele occurs.  Must not be {@code null}.
     * @param variantStart Start position (1-based, inclusive) of the given {@code altAllele}.  Must be > 0.
     * @param variantEnd End position (1-based, inclusive) of the given {@code altAllele}.  Must be > 0.
     * @param exonPositionList List of exon start / stop positions to cross-reference with the given {@code altAllele}.  Must not be {@code null}.
     * @param strand The {@link Strand} on which the {@code altAllele} is located.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return A {@link SimpleInterval} containing the extants of the exon that overlaps the given altAllele.  {@code null} if no overlap occurs.
     */
    public static SimpleInterval getOverlappingExonPositions(final Allele refAllele,
                                                             final Allele altAllele,
                                                             final String contig,
                                                             final int variantStart,
                                                             final int variantEnd,
                                                             final Strand strand,
                                                             final List<? extends htsjdk.samtools.util.Locatable> exonPositionList) {
        Utils.nonNull( refAllele );
        Utils.nonNull( altAllele );
        Utils.nonNull( contig );
        Utils.nonNull( exonPositionList );
        assertValidStrand( strand );

        ParamUtils.isPositive( variantStart, "Genome positions must be > 0." );
        ParamUtils.isPositive( variantEnd, "Genome positions must be > 0." );

        // Get the correct start / end positions:
        int alleleStart = variantStart;
        int alleleEnd   = variantEnd;

        if ( refAllele.length() > refAllele.length() ) {
            final int lengthAdjustment = (refAllele.length() - altAllele.length());
            if ( strand == Strand.NEGATIVE ) {
                alleleStart -= lengthAdjustment;
            }
            else {
                alleleEnd   += lengthAdjustment;
            }
        }

        // Create an interval so we can check for overlap:
        final SimpleInterval variantInterval = new SimpleInterval( contig, alleleStart, alleleEnd );

        int exonStart = -1;
        int exonStop  = -1;

        // Check for overlap:
        for ( final Locatable exon : exonPositionList ) {
            if ( variantInterval.overlaps(exon) ) {
                exonStart = exon.getStart();
                exonStop  = exon.getEnd();
            }
        }

        // Since we set the start/stop together we only need to check one of these:
        if ( exonStart == -1 ) {
            return null;
        }
        else {
            return new SimpleInterval( contig, exonStart, exonStop );
        }
    }

    /**
     * Asserts that the given strand is non-null and is not equal to {@link Strand#NONE}.
     * Throws an Exception if the strand is {@code null} or equal to {@link Strand#NONE}.
     * @param strand The {@link Strand} to validate.
     */
    public static void assertValidStrand( final Strand strand ) {

        Utils.nonNull( strand );

        if ( strand == Strand.NONE ) {
            throw new GATKException("Unable to handle NONE strand.");
        }

    }

    // ========================================================================================

    /**
     * Class representing exceptions that arise when trying to create a coding sequence for a variant:
     */
    public static class TranscriptCodingSequenceException extends GATKException {

        private static final long serialVersionUID = 1L;

        public TranscriptCodingSequenceException( final String msg ) {
            super(msg);
        }

        public TranscriptCodingSequenceException( final String msg, final Throwable throwable ) {
            super(msg, throwable);
        }
    }
}
