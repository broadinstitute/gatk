package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
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
import java.util.stream.Collectors;

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

    private static final SAMSequenceDictionary B37_SEQUENCE_DICTIONARY;

    private static final Map<String, String> B37_To_HG19_CONTIG_NAME_MAP;

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

        B37_SEQUENCE_DICTIONARY = initializeB37SequenceDict();

        B37_To_HG19_CONTIG_NAME_MAP = initializeB37ToHg19ContigNameMap();
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
    public static String createAminoAcidSequence( final String codingSequence ) {
        return createAminoAcidSequence(codingSequence, false);
    }

    /**
     * Creates an amino acid sequence from a given coding sequence.
     * If the coding sequence is not evenly divisible by 3, the remainder bases will not be included in the coding sequence.
     * @param codingSequence The coding sequence from which to create an amino acid sequence.  Must not be {@code null}.
     * @param isFrameshift Whether the given {@code codingSequence} was derived from a frameshift mutation.  In this case, no warning will be issued for incorrect sequence length.
     * @return A {@link String} containing a sequence of single-letter amino acids.
     */
    public static String createAminoAcidSequence(final String codingSequence, final boolean isFrameshift) {

        Utils.nonNull(codingSequence);

        final StringBuilder sb = new StringBuilder();

        // Ensure that we don't have remainder bases:
        int maxIndex = codingSequence.length();
        if ( maxIndex % 3 != 0 ) {
            maxIndex = (int)Math.floor(maxIndex / 3) * 3;
            if ( !isFrameshift ) {
                logger.warn("createAminoAcidSequence given a coding sequence of length not divisible by 3.  Dropping bases from the end: " + (codingSequence.length() % 3));
            }
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
        // TODO: This is generating too long a string for INDELs because of VCF format alleles - see Issue #4407

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
     * Get the position (1-based, inclusive) of the given {@link VariantContext} start relative to the transcript it appears in.
     * The transcript is specified by {@code sortedTranscriptExonList}.
     * @param variant The {@link VariantContext} of which to find the start position in the given transcript (must not be {@code null}).
     * @param exons {@link List} of {@Link Locatable}s representing the exons in the transcript in which the given {@code variant} occurs.
     * @param strand The {@link Strand} on which the {@code variant} occurs.
     * @return The start position (1-based, inclusive) of the given {@code variant} in the transcript in which it appears.
     */
    public static int getTranscriptAlleleStartPosition(final VariantContext variant,
                                                       final List<? extends Locatable> exons,
                                                       final Strand strand) {
        Utils.nonNull(variant);
        Utils.nonNull(exons);
        assertValidStrand(strand);

        // Set up our position variable:
        int position;

        // NOTE: We don't need to worry about UTRs in here - all UTRs occur somewhere in an exon in GENCODE.

        // Filter the elements by whether they come before the variant in the transcript and
        // then sort them by their order in the transcript:
        final List<Locatable> sortedFilteredExons;
        if ( strand == Strand.POSITIVE) {
            sortedFilteredExons = exons.stream()
                    .filter(e -> e.getStart() <= variant.getStart())
                    .sorted(Comparator.comparingInt(Locatable::getStart))
                    .collect(Collectors.toList());

            // We are guaranteed that the variant occurs in the last element of sortedTranscriptElements because of the sorting.
            // Add 1 to position to account for inclusive values:
            position = variant.getStart() - sortedFilteredExons.get(sortedFilteredExons.size()-1).getStart() + 1;
        }
        else {
            sortedFilteredExons = exons.stream()
                    .filter(e -> e.getEnd() >= variant.getStart())
                    .sorted(Comparator.comparingInt(Locatable::getStart).reversed())
                    .collect(Collectors.toList());

            // We are guaranteed that the variant occurs in the last element of sortedTranscriptElements because of the sorting.
            // Add 1 to position to account for inclusive values:
            position = sortedFilteredExons.get(sortedFilteredExons.size()-1).getEnd() - variant.getStart() + 1;
        }

        // Add up the lengths of all exons before the last one:
        for ( int i = 0; i < sortedFilteredExons.size() - 1; ++i ) {
            // Add 1 to position to account for inclusive values:
            position += sortedFilteredExons.get(i).getEnd() - sortedFilteredExons.get(i).getStart() + 1;
        }

        return position;
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
     * Returns if the given {@link SAMSequenceDictionary} is for the B37 Human Genome Reference.
     * @param sequenceDictionary The {@link SAMSequenceDictionary} to check for B37 Reference.
     * @return {@code true} if {@code sequenceDictionary} is for the B37 Human Genome Reference; {@code false} otherwise.
     */
    public static boolean isSequenceDictionaryUsingB37Reference(final SAMSequenceDictionary sequenceDictionary) {
        // Check to make sure all our sequences are accounted for in the given dictionary.

        if ( sequenceDictionary == null ) {
            return false;
        }

        for ( final SAMSequenceRecord b37SequenceRecord : B37_SEQUENCE_DICTIONARY.getSequences() ) {
            // Now we check the Name, Length, and MD5Sum (if present) of all records:

            final SAMSequenceRecord inputSequenceRecord = sequenceDictionary.getSequence(b37SequenceRecord.getSequenceName());
            if ( inputSequenceRecord == null ) {
                return false;
            }

            if ( inputSequenceRecord.getSequenceLength() != b37SequenceRecord.getSequenceLength() ) {
                return false;
            }

            if ( (inputSequenceRecord.getMd5() != null) && (!inputSequenceRecord.getMd5().equals(b37SequenceRecord.getMd5())) ) {
                return false;
            }
        }

        return true;
    }

    /**
     * Converts a given B37 style contig name to the equivalent in hg19.
     * @param b37Contig The contig name from the B37 Human Genome reference to convert to the equivalent contig from the HG19 Human Genome reference.
     * @return The HG19 equivalent of the given B37 contig name, if such an equivalent name exists.  If no equivalent name exists, returns the given {@code b37Contig}.
     */
    public static String convertB37ContigToHg19Contig( final String b37Contig ) {
        if ( B37_To_HG19_CONTIG_NAME_MAP.containsKey(b37Contig) ) {
            return B37_To_HG19_CONTIG_NAME_MAP.get(b37Contig);
        }
        // Couldn't find it in our map, so we just return the original:
        return b37Contig;
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

    /**
     * Get the strand-corrected (reverse complemented) {@link Allele} for the given {@link Allele} and {@link Strand}.
     * @param allele The {@link Allele} to correct for strandedness.
     * @param strand The {@link Strand} on which the given {@code allele} lies.  Must be valid as per {@link #assertValidStrand(Strand)}
     * @return The {@link Allele} with sequence corrected for strand.
     */
    public static Allele getStrandCorrectedAllele(final Allele allele, final Strand strand) {
        assertValidStrand(strand);

        if ( strand == Strand.POSITIVE ) {
            return Allele.create(allele, false);
        }
        else {
            return Allele.create(ReadUtils.getBasesReverseComplement(allele.getBases()), false);
        }
    }

    /**
     * @return An initialized {@link Map} of {@link String} to {@link String} with keys of contig names from the B37 Human Genome Reference and values of the corresponding contig name in the HG19 Human Genome reference.
     */
    private static final Map<String, String> initializeB37ToHg19ContigNameMap() {
        final Map<String, String> b37ToHg19ContigNameMap = new HashMap<>();

        b37ToHg19ContigNameMap.put( "1", "chr1" );
        b37ToHg19ContigNameMap.put( "2", "chr2" );
        b37ToHg19ContigNameMap.put( "3", "chr3" );
        b37ToHg19ContigNameMap.put( "4", "chr4" );
        b37ToHg19ContigNameMap.put( "5", "chr5" );
        b37ToHg19ContigNameMap.put( "6", "chr6" );
        b37ToHg19ContigNameMap.put( "7", "chr7" );
        b37ToHg19ContigNameMap.put( "8", "chr8" );
        b37ToHg19ContigNameMap.put( "9", "chr9" );
        b37ToHg19ContigNameMap.put( "10", "chr10" );
        b37ToHg19ContigNameMap.put( "11", "chr11" );
        b37ToHg19ContigNameMap.put( "12", "chr12" );
        b37ToHg19ContigNameMap.put( "13", "chr13" );
        b37ToHg19ContigNameMap.put( "14", "chr14" );
        b37ToHg19ContigNameMap.put( "15", "chr15" );
        b37ToHg19ContigNameMap.put( "16", "chr16" );
        b37ToHg19ContigNameMap.put( "17", "chr17" );
        b37ToHg19ContigNameMap.put( "18", "chr18" );
        b37ToHg19ContigNameMap.put( "19", "chr19" );
        b37ToHg19ContigNameMap.put( "20", "chr20" );
        b37ToHg19ContigNameMap.put( "21", "chr21" );
        b37ToHg19ContigNameMap.put( "22", "chr22" );
        b37ToHg19ContigNameMap.put( "X", "chrX" );
        b37ToHg19ContigNameMap.put( "Y", "chrY" );
        b37ToHg19ContigNameMap.put( "MT", "chrM" );

        // These are included for completeness, but are not necessary as we default to the input contig name if the
        // contig is not present in this map.  (Because they are the same between both.)
        b37ToHg19ContigNameMap.put( "GL000202.1", "GL000202.1" );
        b37ToHg19ContigNameMap.put( "GL000244.1", "GL000244.1" );
        b37ToHg19ContigNameMap.put( "GL000235.1", "GL000235.1" );
        b37ToHg19ContigNameMap.put( "GL000238.1", "GL000238.1" );
        b37ToHg19ContigNameMap.put( "GL000226.1", "GL000226.1" );
        b37ToHg19ContigNameMap.put( "GL000218.1", "GL000218.1" );
        b37ToHg19ContigNameMap.put( "GL000249.1", "GL000249.1" );
        b37ToHg19ContigNameMap.put( "GL000242.1", "GL000242.1" );
        b37ToHg19ContigNameMap.put( "GL000221.1", "GL000221.1" );
        b37ToHg19ContigNameMap.put( "GL000192.1", "GL000192.1" );
        b37ToHg19ContigNameMap.put( "GL000223.1", "GL000223.1" );
        b37ToHg19ContigNameMap.put( "GL000232.1", "GL000232.1" );
        b37ToHg19ContigNameMap.put( "GL000206.1", "GL000206.1" );
        b37ToHg19ContigNameMap.put( "GL000240.1", "GL000240.1" );
        b37ToHg19ContigNameMap.put( "GL000214.1", "GL000214.1" );
        b37ToHg19ContigNameMap.put( "GL000212.1", "GL000212.1" );
        b37ToHg19ContigNameMap.put( "GL000199.1", "GL000199.1" );
        b37ToHg19ContigNameMap.put( "GL000248.1", "GL000248.1" );
        b37ToHg19ContigNameMap.put( "GL000195.1", "GL000195.1" );
        b37ToHg19ContigNameMap.put( "GL000215.1", "GL000215.1" );
        b37ToHg19ContigNameMap.put( "GL000225.1", "GL000225.1" );
        b37ToHg19ContigNameMap.put( "GL000216.1", "GL000216.1" );
        b37ToHg19ContigNameMap.put( "GL000194.1", "GL000194.1" );
        b37ToHg19ContigNameMap.put( "GL000217.1", "GL000217.1" );
        b37ToHg19ContigNameMap.put( "GL000197.1", "GL000197.1" );
        b37ToHg19ContigNameMap.put( "GL000222.1", "GL000222.1" );
        b37ToHg19ContigNameMap.put( "GL000200.1", "GL000200.1" );
        b37ToHg19ContigNameMap.put( "GL000211.1", "GL000211.1" );
        b37ToHg19ContigNameMap.put( "GL000247.1", "GL000247.1" );
        b37ToHg19ContigNameMap.put( "GL000233.1", "GL000233.1" );
        b37ToHg19ContigNameMap.put( "GL000210.1", "GL000210.1" );
        b37ToHg19ContigNameMap.put( "GL000198.1", "GL000198.1" );
        b37ToHg19ContigNameMap.put( "GL000245.1", "GL000245.1" );
        b37ToHg19ContigNameMap.put( "GL000234.1", "GL000234.1" );
        b37ToHg19ContigNameMap.put( "GL000203.1", "GL000203.1" );
        b37ToHg19ContigNameMap.put( "GL000239.1", "GL000239.1" );
        b37ToHg19ContigNameMap.put( "GL000213.1", "GL000213.1" );
        b37ToHg19ContigNameMap.put( "GL000227.1", "GL000227.1" );
        b37ToHg19ContigNameMap.put( "GL000208.1", "GL000208.1" );
        b37ToHg19ContigNameMap.put( "GL000230.1", "GL000230.1" );
        b37ToHg19ContigNameMap.put( "GL000231.1", "GL000231.1" );
        b37ToHg19ContigNameMap.put( "GL000228.1", "GL000228.1" );
        b37ToHg19ContigNameMap.put( "GL000243.1", "GL000243.1" );
        b37ToHg19ContigNameMap.put( "GL000229.1", "GL000229.1" );
        b37ToHg19ContigNameMap.put( "GL000205.1", "GL000205.1" );
        b37ToHg19ContigNameMap.put( "GL000224.1", "GL000224.1" );
        b37ToHg19ContigNameMap.put( "GL000191.1", "GL000191.1" );
        b37ToHg19ContigNameMap.put( "GL000196.1", "GL000196.1" );
        b37ToHg19ContigNameMap.put( "GL000193.1", "GL000193.1" );
        b37ToHg19ContigNameMap.put( "GL000201.1", "GL000201.1" );
        b37ToHg19ContigNameMap.put( "GL000237.1", "GL000237.1" );
        b37ToHg19ContigNameMap.put( "GL000246.1", "GL000246.1" );
        b37ToHg19ContigNameMap.put( "GL000241.1", "GL000241.1" );
        b37ToHg19ContigNameMap.put( "GL000204.1", "GL000204.1" );
        b37ToHg19ContigNameMap.put( "GL000207.1", "GL000207.1" );
        b37ToHg19ContigNameMap.put( "GL000209.1", "GL000209.1" );
        b37ToHg19ContigNameMap.put( "GL000219.1", "GL000219.1" );
        b37ToHg19ContigNameMap.put( "GL000220.1", "GL000220.1" );
        b37ToHg19ContigNameMap.put( "GL000236.1", "GL000236.1" );

        return b37ToHg19ContigNameMap;
    }

    /**
     * @return An initialized {@link SAMSequenceDictionary} containing the {@link SAMSequenceRecord}s from the B37 Human Genome Reference.
     */
    private static final SAMSequenceDictionary initializeB37SequenceDict() {
        final SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary();

        final SAMSequenceRecord seq1 = new SAMSequenceRecord("1", 249250621);
        seq1.setAssembly("GRCh37");
        seq1.setMd5("1b22b98cdeb4a9304cb5d48026a85128");
        seq1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq1);
        
        final SAMSequenceRecord seq2 = new SAMSequenceRecord( "2", 243199373 );
        seq2.setAssembly("GRCh37");
        seq2.setMd5("a0d9851da00400dec1098a9255ac712e");
        seq2.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq2);
        
        final SAMSequenceRecord seq3 = new SAMSequenceRecord( "3", 198022430 );
        seq3.setAssembly("GRCh37");
        seq3.setMd5("fdfd811849cc2fadebc929bb925902e5");
        seq3.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq3);
        
        final SAMSequenceRecord seq4 = new SAMSequenceRecord( "4", 191154276 );
        seq4.setAssembly("GRCh37");
        seq4.setMd5("23dccd106897542ad87d2765d28a19a1");
        seq4.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq4);
        
        final SAMSequenceRecord seq5 = new SAMSequenceRecord( "5", 180915260 );
        seq5.setAssembly("GRCh37");
        seq5.setMd5("0740173db9ffd264d728f32784845cd7");
        seq5.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq5);
        
        final SAMSequenceRecord seq6 = new SAMSequenceRecord( "6", 171115067 );
        seq6.setAssembly("GRCh37");
        seq6.setMd5("1d3a93a248d92a729ee764823acbbc6b");
        seq6.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq6);
        
        final SAMSequenceRecord seq7 = new SAMSequenceRecord( "7", 159138663 );
        seq7.setAssembly("GRCh37");
        seq7.setMd5("618366e953d6aaad97dbe4777c29375e");
        seq7.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq7);
        
        final SAMSequenceRecord seq8 = new SAMSequenceRecord( "8", 146364022 );
        seq8.setAssembly("GRCh37");
        seq8.setMd5("96f514a9929e410c6651697bded59aec");
        seq8.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq8);
        
        final SAMSequenceRecord seq9 = new SAMSequenceRecord( "9", 141213431 );
        seq9.setAssembly("GRCh37");
        seq9.setMd5("3e273117f15e0a400f01055d9f393768");
        seq9.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq9);
        
        final SAMSequenceRecord seq10 = new SAMSequenceRecord( "10", 135534747 );
        seq10.setAssembly("GRCh37");
        seq10.setMd5("988c28e000e84c26d552359af1ea2e1d");
        seq10.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq10);
        
        final SAMSequenceRecord seq11 = new SAMSequenceRecord( "11", 135006516 );
        seq11.setAssembly("GRCh37");
        seq11.setMd5("98c59049a2df285c76ffb1c6db8f8b96");
        seq11.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq11);
        
        final SAMSequenceRecord seq12 = new SAMSequenceRecord( "12", 133851895 );
        seq12.setAssembly("GRCh37");
        seq12.setMd5("51851ac0e1a115847ad36449b0015864");
        seq12.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq12);
        
        final SAMSequenceRecord seq13 = new SAMSequenceRecord( "13", 115169878 );
        seq13.setAssembly("GRCh37");
        seq13.setMd5("283f8d7892baa81b510a015719ca7b0b");
        seq13.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq13);
        
        final SAMSequenceRecord seq14 = new SAMSequenceRecord( "14", 107349540 );
        seq14.setAssembly("GRCh37");
        seq14.setMd5("98f3cae32b2a2e9524bc19813927542e");
        seq14.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq14);
        
        final SAMSequenceRecord seq15 = new SAMSequenceRecord( "15", 102531392 );
        seq15.setAssembly("GRCh37");
        seq15.setMd5("e5645a794a8238215b2cd77acb95a078");
        seq15.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq15);
        
        final SAMSequenceRecord seq16 = new SAMSequenceRecord( "16", 90354753 );
        seq16.setAssembly("GRCh37");
        seq16.setMd5("fc9b1a7b42b97a864f56b348b06095e6");
        seq16.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq16);
        
        final SAMSequenceRecord seq17 = new SAMSequenceRecord( "17", 81195210 );
        seq17.setAssembly("GRCh37");
        seq17.setMd5("351f64d4f4f9ddd45b35336ad97aa6de");
        seq17.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq17);
        
        final SAMSequenceRecord seq18 = new SAMSequenceRecord( "18", 78077248 );
        seq18.setAssembly("GRCh37");
        seq18.setMd5("b15d4b2d29dde9d3e4f93d1d0f2cbc9c");
        seq18.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq18);
        
        final SAMSequenceRecord seq19 = new SAMSequenceRecord( "19", 59128983 );
        seq19.setAssembly("GRCh37");
        seq19.setMd5("1aacd71f30db8e561810913e0b72636d");
        seq19.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq19);
        
        final SAMSequenceRecord seq20 = new SAMSequenceRecord( "20", 63025520 );
        seq20.setAssembly("GRCh37");
        seq20.setMd5("0dec9660ec1efaaf33281c0d5ea2560f");
        seq20.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq20);
        
        final SAMSequenceRecord seq21 = new SAMSequenceRecord( "21", 48129895 );
        seq21.setAssembly("GRCh37");
        seq21.setMd5("2979a6085bfe28e3ad6f552f361ed74d");
        seq21.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq21);
        
        final SAMSequenceRecord seq22 = new SAMSequenceRecord( "22", 51304566 );
        seq22.setAssembly("GRCh37");
        seq22.setMd5("a718acaa6135fdca8357d5bfe94211dd");
        seq22.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seq22);
        
        final SAMSequenceRecord seqX = new SAMSequenceRecord( "X", 155270560 );
        seqX.setAssembly("GRCh37");
        seqX.setMd5("7e0e2e580297b7764e31dbc80c2540dd");
        seqX.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqX);
        
        final SAMSequenceRecord seqY = new SAMSequenceRecord( "Y", 59373566 );
        seqY.setAssembly("GRCh37");
        seqY.setMd5("1fa3474750af0948bdf97d5a0ee52e51");
        seqY.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqY);
        
        final SAMSequenceRecord seqMT = new SAMSequenceRecord( "MT", 16569 );
        seqMT.setAssembly("GRCh37");
        seqMT.setMd5("c68f52674c9fb33aef52dcf399755519");
        seqMT.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqMT);
        
        final SAMSequenceRecord seqGL000207_1 = new SAMSequenceRecord( "GL000207.1", 4262 );
        seqGL000207_1.setAssembly("GRCh37");
        seqGL000207_1.setMd5("f3814841f1939d3ca19072d9e89f3fd7");
        seqGL000207_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000207_1);
        
        final SAMSequenceRecord seqGL000226_1 = new SAMSequenceRecord( "GL000226.1", 15008 );
        seqGL000226_1.setAssembly("GRCh37");
        seqGL000226_1.setMd5("1c1b2cd1fccbc0a99b6a447fa24d1504");
        seqGL000226_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000226_1);
        
        final SAMSequenceRecord seqGL000229_1 = new SAMSequenceRecord( "GL000229.1", 19913 );
        seqGL000229_1.setAssembly("GRCh37");
        seqGL000229_1.setMd5("d0f40ec87de311d8e715b52e4c7062e1");
        seqGL000229_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000229_1);
        
        final SAMSequenceRecord seqGL000231_1 = new SAMSequenceRecord( "GL000231.1", 27386 );
        seqGL000231_1.setAssembly("GRCh37");
        seqGL000231_1.setMd5("ba8882ce3a1efa2080e5d29b956568a4");
        seqGL000231_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000231_1);
        
        final SAMSequenceRecord seqGL000210_1 = new SAMSequenceRecord( "GL000210.1", 27682 );
        seqGL000210_1.setAssembly("GRCh37");
        seqGL000210_1.setMd5("851106a74238044126131ce2a8e5847c");
        seqGL000210_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000210_1);
        
        final SAMSequenceRecord seqGL000239_1 = new SAMSequenceRecord( "GL000239.1", 33824 );
        seqGL000239_1.setAssembly("GRCh37");
        seqGL000239_1.setMd5("99795f15702caec4fa1c4e15f8a29c07");
        seqGL000239_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000239_1);
        
        final SAMSequenceRecord seqGL000235_1 = new SAMSequenceRecord( "GL000235.1", 34474 );
        seqGL000235_1.setAssembly("GRCh37");
        seqGL000235_1.setMd5("118a25ca210cfbcdfb6c2ebb249f9680");
        seqGL000235_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000235_1);
        
        final SAMSequenceRecord seqGL000201_1 = new SAMSequenceRecord( "GL000201.1", 36148 );
        seqGL000201_1.setAssembly("GRCh37");
        seqGL000201_1.setMd5("dfb7e7ec60ffdcb85cb359ea28454ee9");
        seqGL000201_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000201_1);
        
        final SAMSequenceRecord seqGL000247_1 = new SAMSequenceRecord( "GL000247.1", 36422 );
        seqGL000247_1.setAssembly("GRCh37");
        seqGL000247_1.setMd5("7de00226bb7df1c57276ca6baabafd15");
        seqGL000247_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000247_1);
        
        final SAMSequenceRecord seqGL000245_1 = new SAMSequenceRecord( "GL000245.1", 36651 );
        seqGL000245_1.setAssembly("GRCh37");
        seqGL000245_1.setMd5("89bc61960f37d94abf0df2d481ada0ec");
        seqGL000245_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000245_1);
        
        final SAMSequenceRecord seqGL000197_1 = new SAMSequenceRecord( "GL000197.1", 37175 );
        seqGL000197_1.setAssembly("GRCh37");
        seqGL000197_1.setMd5("6f5efdd36643a9b8c8ccad6f2f1edc7b");
        seqGL000197_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000197_1);
        
        final SAMSequenceRecord seqGL000203_1 = new SAMSequenceRecord( "GL000203.1", 37498 );
        seqGL000203_1.setAssembly("GRCh37");
        seqGL000203_1.setMd5("96358c325fe0e70bee73436e8bb14dbd");
        seqGL000203_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000203_1);
        
        final SAMSequenceRecord seqGL000246_1 = new SAMSequenceRecord( "GL000246.1", 38154 );
        seqGL000246_1.setAssembly("GRCh37");
        seqGL000246_1.setMd5("e4afcd31912af9d9c2546acf1cb23af2");
        seqGL000246_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000246_1);
        
        final SAMSequenceRecord seqGL000249_1 = new SAMSequenceRecord( "GL000249.1", 38502 );
        seqGL000249_1.setAssembly("GRCh37");
        seqGL000249_1.setMd5("1d78abec37c15fe29a275eb08d5af236");
        seqGL000249_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000249_1);
        
        final SAMSequenceRecord seqGL000196_1 = new SAMSequenceRecord( "GL000196.1", 38914 );
        seqGL000196_1.setAssembly("GRCh37");
        seqGL000196_1.setMd5("d92206d1bb4c3b4019c43c0875c06dc0");
        seqGL000196_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000196_1);
        
        final SAMSequenceRecord seqGL000248_1 = new SAMSequenceRecord( "GL000248.1", 39786 );
        seqGL000248_1.setAssembly("GRCh37");
        seqGL000248_1.setMd5("5a8e43bec9be36c7b49c84d585107776");
        seqGL000248_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000248_1);
        
        final SAMSequenceRecord seqGL000244_1 = new SAMSequenceRecord( "GL000244.1", 39929 );
        seqGL000244_1.setAssembly("GRCh37");
        seqGL000244_1.setMd5("0996b4475f353ca98bacb756ac479140");
        seqGL000244_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000244_1);
        
        final SAMSequenceRecord seqGL000238_1 = new SAMSequenceRecord( "GL000238.1", 39939 );
        seqGL000238_1.setAssembly("GRCh37");
        seqGL000238_1.setMd5("131b1efc3270cc838686b54e7c34b17b");
        seqGL000238_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000238_1);
        
        final SAMSequenceRecord seqGL000202_1 = new SAMSequenceRecord( "GL000202.1", 40103 );
        seqGL000202_1.setAssembly("GRCh37");
        seqGL000202_1.setMd5("06cbf126247d89664a4faebad130fe9c");
        seqGL000202_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000202_1);
        
        final SAMSequenceRecord seqGL000234_1 = new SAMSequenceRecord( "GL000234.1", 40531 );
        seqGL000234_1.setAssembly("GRCh37");
        seqGL000234_1.setMd5("93f998536b61a56fd0ff47322a911d4b");
        seqGL000234_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000234_1);
        
        final SAMSequenceRecord seqGL000232_1 = new SAMSequenceRecord( "GL000232.1", 40652 );
        seqGL000232_1.setAssembly("GRCh37");
        seqGL000232_1.setMd5("3e06b6741061ad93a8587531307057d8");
        seqGL000232_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000232_1);
        
        final SAMSequenceRecord seqGL000206_1 = new SAMSequenceRecord( "GL000206.1", 41001 );
        seqGL000206_1.setAssembly("GRCh37");
        seqGL000206_1.setMd5("43f69e423533e948bfae5ce1d45bd3f1");
        seqGL000206_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000206_1);
        
        final SAMSequenceRecord seqGL000240_1 = new SAMSequenceRecord( "GL000240.1", 41933 );
        seqGL000240_1.setAssembly("GRCh37");
        seqGL000240_1.setMd5("445a86173da9f237d7bcf41c6cb8cc62");
        seqGL000240_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000240_1);
        
        final SAMSequenceRecord seqGL000236_1 = new SAMSequenceRecord( "GL000236.1", 41934 );
        seqGL000236_1.setAssembly("GRCh37");
        seqGL000236_1.setMd5("fdcd739913efa1fdc64b6c0cd7016779");
        seqGL000236_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000236_1);
        
        final SAMSequenceRecord seqGL000241_1 = new SAMSequenceRecord( "GL000241.1", 42152 );
        seqGL000241_1.setAssembly("GRCh37");
        seqGL000241_1.setMd5("ef4258cdc5a45c206cea8fc3e1d858cf");
        seqGL000241_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000241_1);
        
        final SAMSequenceRecord seqGL000243_1 = new SAMSequenceRecord( "GL000243.1", 43341 );
        seqGL000243_1.setAssembly("GRCh37");
        seqGL000243_1.setMd5("cc34279a7e353136741c9fce79bc4396");
        seqGL000243_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000243_1);
        
        final SAMSequenceRecord seqGL000242_1 = new SAMSequenceRecord( "GL000242.1", 43523 );
        seqGL000242_1.setAssembly("GRCh37");
        seqGL000242_1.setMd5("2f8694fc47576bc81b5fe9e7de0ba49e");
        seqGL000242_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000242_1);
        
        final SAMSequenceRecord seqGL000230_1 = new SAMSequenceRecord( "GL000230.1", 43691 );
        seqGL000230_1.setAssembly("GRCh37");
        seqGL000230_1.setMd5("b4eb71ee878d3706246b7c1dbef69299");
        seqGL000230_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000230_1);
        
        final SAMSequenceRecord seqGL000237_1 = new SAMSequenceRecord( "GL000237.1", 45867 );
        seqGL000237_1.setAssembly("GRCh37");
        seqGL000237_1.setMd5("e0c82e7751df73f4f6d0ed30cdc853c0");
        seqGL000237_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000237_1);
        
        final SAMSequenceRecord seqGL000233_1 = new SAMSequenceRecord( "GL000233.1", 45941 );
        seqGL000233_1.setAssembly("GRCh37");
        seqGL000233_1.setMd5("7fed60298a8d62ff808b74b6ce820001");
        seqGL000233_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000233_1);
        
        final SAMSequenceRecord seqGL000204_1 = new SAMSequenceRecord( "GL000204.1", 81310 );
        seqGL000204_1.setAssembly("GRCh37");
        seqGL000204_1.setMd5("efc49c871536fa8d79cb0a06fa739722");
        seqGL000204_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000204_1);
        
        final SAMSequenceRecord seqGL000198_1 = new SAMSequenceRecord( "GL000198.1", 90085 );
        seqGL000198_1.setAssembly("GRCh37");
        seqGL000198_1.setMd5("868e7784040da90d900d2d1b667a1383");
        seqGL000198_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000198_1);
        
        final SAMSequenceRecord seqGL000208_1 = new SAMSequenceRecord( "GL000208.1", 92689 );
        seqGL000208_1.setAssembly("GRCh37");
        seqGL000208_1.setMd5("aa81be49bf3fe63a79bdc6a6f279abf6");
        seqGL000208_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000208_1);
        
        final SAMSequenceRecord seqGL000191_1 = new SAMSequenceRecord( "GL000191.1", 106433 );
        seqGL000191_1.setAssembly("GRCh37");
        seqGL000191_1.setMd5("d75b436f50a8214ee9c2a51d30b2c2cc");
        seqGL000191_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000191_1);
        
        final SAMSequenceRecord seqGL000227_1 = new SAMSequenceRecord( "GL000227.1", 128374 );
        seqGL000227_1.setAssembly("GRCh37");
        seqGL000227_1.setMd5("a4aead23f8053f2655e468bcc6ecdceb");
        seqGL000227_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000227_1);
        
        final SAMSequenceRecord seqGL000228_1 = new SAMSequenceRecord( "GL000228.1", 129120 );
        seqGL000228_1.setAssembly("GRCh37");
        seqGL000228_1.setMd5("c5a17c97e2c1a0b6a9cc5a6b064b714f");
        seqGL000228_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000228_1);
        
        final SAMSequenceRecord seqGL000214_1 = new SAMSequenceRecord( "GL000214.1", 137718 );
        seqGL000214_1.setAssembly("GRCh37");
        seqGL000214_1.setMd5("46c2032c37f2ed899eb41c0473319a69");
        seqGL000214_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000214_1);
        
        final SAMSequenceRecord seqGL000221_1 = new SAMSequenceRecord( "GL000221.1", 155397 );
        seqGL000221_1.setAssembly("GRCh37");
        seqGL000221_1.setMd5("3238fb74ea87ae857f9c7508d315babb");
        seqGL000221_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000221_1);
        
        final SAMSequenceRecord seqGL000209_1 = new SAMSequenceRecord( "GL000209.1", 159169 );
        seqGL000209_1.setAssembly("GRCh37");
        seqGL000209_1.setMd5("f40598e2a5a6b26e84a3775e0d1e2c81");
        seqGL000209_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000209_1);
        
        final SAMSequenceRecord seqGL000218_1 = new SAMSequenceRecord( "GL000218.1", 161147 );
        seqGL000218_1.setAssembly("GRCh37");
        seqGL000218_1.setMd5("1d708b54644c26c7e01c2dad5426d38c");
        seqGL000218_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000218_1);
        
        final SAMSequenceRecord seqGL000220_1 = new SAMSequenceRecord( "GL000220.1", 161802 );
        seqGL000220_1.setAssembly("GRCh37");
        seqGL000220_1.setMd5("fc35de963c57bf7648429e6454f1c9db");
        seqGL000220_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000220_1);
        
        final SAMSequenceRecord seqGL000213_1 = new SAMSequenceRecord( "GL000213.1", 164239 );
        seqGL000213_1.setAssembly("GRCh37");
        seqGL000213_1.setMd5("9d424fdcc98866650b58f004080a992a");
        seqGL000213_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000213_1);
        
        final SAMSequenceRecord seqGL000211_1 = new SAMSequenceRecord( "GL000211.1", 166566 );
        seqGL000211_1.setAssembly("GRCh37");
        seqGL000211_1.setMd5("7daaa45c66b288847b9b32b964e623d3");
        seqGL000211_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000211_1);
        
        final SAMSequenceRecord seqGL000199_1 = new SAMSequenceRecord( "GL000199.1", 169874 );
        seqGL000199_1.setAssembly("GRCh37");
        seqGL000199_1.setMd5("569af3b73522fab4b40995ae4944e78e");
        seqGL000199_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000199_1);
        
        final SAMSequenceRecord seqGL000217_1 = new SAMSequenceRecord( "GL000217.1", 172149 );
        seqGL000217_1.setAssembly("GRCh37");
        seqGL000217_1.setMd5("6d243e18dea1945fb7f2517615b8f52e");
        seqGL000217_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000217_1);
        
        final SAMSequenceRecord seqGL000216_1 = new SAMSequenceRecord( "GL000216.1", 172294 );
        seqGL000216_1.setAssembly("GRCh37");
        seqGL000216_1.setMd5("642a232d91c486ac339263820aef7fe0");
        seqGL000216_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000216_1);
        
        final SAMSequenceRecord seqGL000215_1 = new SAMSequenceRecord( "GL000215.1", 172545 );
        seqGL000215_1.setAssembly("GRCh37");
        seqGL000215_1.setMd5("5eb3b418480ae67a997957c909375a73");
        seqGL000215_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000215_1);
        
        final SAMSequenceRecord seqGL000205_1 = new SAMSequenceRecord( "GL000205.1", 174588 );
        seqGL000205_1.setAssembly("GRCh37");
        seqGL000205_1.setMd5("d22441398d99caf673e9afb9a1908ec5");
        seqGL000205_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000205_1);
        
        final SAMSequenceRecord seqGL000219_1 = new SAMSequenceRecord( "GL000219.1", 179198 );
        seqGL000219_1.setAssembly("GRCh37");
        seqGL000219_1.setMd5("f977edd13bac459cb2ed4a5457dba1b3");
        seqGL000219_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000219_1);
        
        final SAMSequenceRecord seqGL000224_1 = new SAMSequenceRecord( "GL000224.1", 179693 );
        seqGL000224_1.setAssembly("GRCh37");
        seqGL000224_1.setMd5("d5b2fc04f6b41b212a4198a07f450e20");
        seqGL000224_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000224_1);
        
        final SAMSequenceRecord seqGL000223_1 = new SAMSequenceRecord( "GL000223.1", 180455 );
        seqGL000223_1.setAssembly("GRCh37");
        seqGL000223_1.setMd5("399dfa03bf32022ab52a846f7ca35b30");
        seqGL000223_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000223_1);
        
        final SAMSequenceRecord seqGL000195_1 = new SAMSequenceRecord( "GL000195.1", 182896 );
        seqGL000195_1.setAssembly("GRCh37");
        seqGL000195_1.setMd5("5d9ec007868d517e73543b005ba48535");
        seqGL000195_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000195_1);
        
        final SAMSequenceRecord seqGL000212_1 = new SAMSequenceRecord( "GL000212.1", 186858 );
        seqGL000212_1.setAssembly("GRCh37");
        seqGL000212_1.setMd5("563531689f3dbd691331fd6c5730a88b");
        seqGL000212_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000212_1);
        
        final SAMSequenceRecord seqGL000222_1 = new SAMSequenceRecord( "GL000222.1", 186861 );
        seqGL000222_1.setAssembly("GRCh37");
        seqGL000222_1.setMd5("6fe9abac455169f50470f5a6b01d0f59");
        seqGL000222_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000222_1);
        
        final SAMSequenceRecord seqGL000200_1 = new SAMSequenceRecord( "GL000200.1", 187035 );
        seqGL000200_1.setAssembly("GRCh37");
        seqGL000200_1.setMd5("75e4c8d17cd4addf3917d1703cacaf25");
        seqGL000200_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000200_1);
        
        final SAMSequenceRecord seqGL000193_1 = new SAMSequenceRecord( "GL000193.1", 189789 );
        seqGL000193_1.setAssembly("GRCh37");
        seqGL000193_1.setMd5("dbb6e8ece0b5de29da56601613007c2a");
        seqGL000193_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000193_1);
        
        final SAMSequenceRecord seqGL000194_1 = new SAMSequenceRecord( "GL000194.1", 191469 );
        seqGL000194_1.setAssembly("GRCh37");
        seqGL000194_1.setMd5("6ac8f815bf8e845bb3031b73f812c012");
        seqGL000194_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000194_1);
        
        final SAMSequenceRecord seqGL000225_1 = new SAMSequenceRecord( "GL000225.1", 211173 );
        seqGL000225_1.setAssembly("GRCh37");
        seqGL000225_1.setMd5("63945c3e6962f28ffd469719a747e73c");
        seqGL000225_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000225_1);
        
        final SAMSequenceRecord seqGL000192_1 = new SAMSequenceRecord( "GL000192.1", 547496 );
        seqGL000192_1.setAssembly("GRCh37");
        seqGL000192_1.setMd5("325ba9e808f669dfeee210fdd7b470ac");
        seqGL000192_1.setSpecies("Homo Sapiens");
        sequenceDictionary.addSequence(seqGL000192_1);
        
        final SAMSequenceRecord seqNC_007605 = new SAMSequenceRecord( "NC_007605", 171823 );
        seqNC_007605.setAssembly("NC_007605.1");
        seqNC_007605.setMd5("6743bd63b3ff2b5b8985d8933c53290a");
        seqNC_007605.setSpecies("Epstein-Barr virus");
        sequenceDictionary.addSequence(seqNC_007605);

        return sequenceDictionary;
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
