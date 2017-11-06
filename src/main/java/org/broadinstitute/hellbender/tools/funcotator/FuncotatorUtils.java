/*
* Copyright 2012-2016 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.hellbender.tools.funcotator;

/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class FuncotatorUtils {

    private final static Logger logger = Logger.getLogger(FuncotatorUtils.class);

    /**
     * PRIVATE CONSTRUCTOR
     * DO NOT INSTANTIATE THIS CLASS!
     */
    private FuncotatorUtils() {}

    private static final HashMap<String, AminoAcid> tableByCodon = new HashMap<>(AminoAcid.values().length);
    private static final HashMap<String, AminoAcid> tableByCode = new HashMap<>(AminoAcid.values().length);
    private static final HashMap<String, AminoAcid> tableByLetter = new HashMap<>(AminoAcid.values().length);

    /**
     * Initialize our hashmaps of lookup tables:
     */
    static {
        for ( final AminoAcid acid : AminoAcid.values() ) {
            tableByCode.put(acid.getCode(),acid);
            tableByLetter.put(acid.getLetter(), acid);
            for ( final String codon : acid.getCodons() ) {
                tableByCodon.put(codon,acid);
            }
        }
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
     * Creates the string representation of the codon change for the given {@link SequenceComparison}.
     * Requires that the given {@code seqComp} has the following fields defined with values that are not {@code null}:
     *     alignedCodingSequenceAlleleStart
     *     alignedReferenceAlleleStop
     *     alignedCodingSequenceReferenceAllele
     *     alternateAllele
     *     alignedCodingSequenceAlternateAllele
     *     codingSequenceAlleleStart
     * @param seqComp {@link SequenceComparison} representing the alternate and reference alleles for a DNA sequence.  Must not be {@code null}.
     * @return A {@link String} representing the codon change for the given {@link SequenceComparison}.
     */
    public static String getCodonChangeString(final SequenceComparison seqComp) {

        Utils.nonNull(seqComp);
        Utils.nonNull(seqComp.getAlignedCodingSequenceAlleleStart());
        Utils.nonNull(seqComp.getAlignedReferenceAlleleStop());
        Utils.nonNull(seqComp.getAlignedCodingSequenceReferenceAllele());
        Utils.nonNull(seqComp.getAlternateAllele());
        Utils.nonNull(seqComp.getAlignedCodingSequenceAlternateAllele());
        Utils.nonNull(seqComp.getCodingSequenceAlleleStart());

        // Used for insertions:
        Utils.nonNull(seqComp.getReferenceCodingSequence());
        assertValidStrand( seqComp.getStrand() );

        // Make some local names so that things are easier to read:
        final String alignedRefAllele = seqComp.getAlignedCodingSequenceReferenceAllele();
        final String alignedAltAllele = seqComp.getAlignedCodingSequenceAlternateAllele();
        
        final StringBuilder ref = new StringBuilder();
        final StringBuilder alt = new StringBuilder();

        // Handle the ONP case:
        if ( GATKProtectedVariantContextUtils.isOnp(alignedRefAllele, alignedAltAllele) ) {
            return getCodonChangeStringForOnp( alignedRefAllele, alignedAltAllele,
                    seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop() );
        }
        // Handle the insertion case:
        else if ( GATKProtectedVariantContextUtils.isInsertion(alignedRefAllele, alignedAltAllele) ) {

            final String nextRefCodon = getNextReferenceCodon(seqComp.getReferenceCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());

            // Check for frame shift syntax first:
            if ( GATKProtectedVariantContextUtils.isFrameshift(alignedRefAllele, alignedAltAllele) ) {
                return getCodonChangeStringForInsertionFrameshift(
                        seqComp.getCodingSequenceAlleleStart(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(),
                        nextRefCodon, seqComp.getAlignedReferenceAllele() );
            }
            else {

                // If the variant position is in frame, then we have to use the next codon, because it will be the first
                // codon that is affected.
                // We add 1 because of the convention that the variant occurs just after the base specified for insertions.
                // TODO: This is really bad - we are tying our output to a specific input format.  FIX IT.
                if (isPositionInFrame(seqComp.getCodingSequenceAlleleStart() + 1)) {
                    return getCodonChangeStringForInsertionInFrameWithInFrameStartPosition(
                            alignedRefAllele,
                            alignedAltAllele,
                            nextRefCodon,
                            seqComp.getAlignedCodingSequenceAlleleStart(),
                            seqComp.getAlignedReferenceAlleleStop()
                    );
                }
                else {
                    return getCodonChangeStringForInsertionInFrameWithOutOfFrameStartPosition(
                            alignedRefAllele,
                            alignedAltAllele,
                            seqComp.getAlternateAllele().length(),
                            seqComp.getAlignedCodingSequenceAlleleStart(),
                            seqComp.getAlignedReferenceAlleleStop()
                    );
                }
            }
        }
        // Handle the deletion case:
        else {
            return getCodonChangeStringForDeletion(
                    alignedRefAllele,
                    alignedAltAllele,
                    seqComp.getAlignedCodingSequenceAlleleStart(),
                    seqComp.getAlignedReferenceAlleleStop()
            );
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
     * Get the codon change string for a Frameshift Insertion.
     * @param codingSequenceAlleleStart The position (1-based, inclusive) of the start of the coding sequence of the alleles.  Not guaranteed to be in-frame.
     * @param alignedCodingSequenceAlleleStart The position (1-based, inclusive) of the start of the first of the alleles in the coding sequence (prepended with up to 2 bases to be in-frame).
     * @param alignedReferenceAlleleStop The position (1-based, inclusive) of the in-frame end of the Reference allele.
     * @param nextRefCodon The {@link String} of 3 bases in the next codon in the reference sequence.
     * @param alignedRefAllele The {@link String} of bases contained in the Reference allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @return The {@link String} containing the codon change string for the given alleleic information.
     */
    private static String getCodonChangeStringForInsertionFrameshift(final int codingSequenceAlleleStart,
                                                                     final int alignedCodingSequenceAlleleStart,
                                                                     final int alignedReferenceAlleleStop,
                                                                     final String nextRefCodon,
                                                                     final String alignedRefAllele) {
        // If the variant position is in frame, then we have to use the next codon, because it will be the first
        // codon that is affected.
        // We add 1 because of the convention that the variant occurs just after the base specified for insertions.
        // TODO: This is really bad - we are tying our output to a specific input format.  FIX IT.
        if ( isPositionInFrame(codingSequenceAlleleStart + 1) ) {
            return "c.(" + (alignedCodingSequenceAlleleStart + 3) + "-" +
                    (alignedReferenceAlleleStop + 3) + ")" +
                    nextRefCodon.toLowerCase() + "fs";
        }
        else {
            // Out of frame insertions affect the current codon, so the current position should be
            // used.
            return "c.(" + alignedCodingSequenceAlleleStart + "-" +
                    alignedReferenceAlleleStop + ")" +
                    alignedRefAllele.toLowerCase() + "fs";
        }
    }

    /**
     * Get the codon change string for an in-frame Insertion with an in-frame start position.
     * @param alignedRefAllele The {@link String} of bases contained in the Reference allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param alignedAltAllele The {@link String} of bases contained in the Alternate allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param nextRefCodon The {@link String} of 3 bases in the next codon in the reference sequence.
     * @param alignedCodingSequenceAlleleStart The position (1-based, inclusive) of the start of the first of the alleles in the coding sequence (prepended with up to 2 bases to be in-frame).
     * @param alignedReferenceAlleleStop The position (1-based, inclusive) of the in-frame end of the Reference allele.
     * @return The {@link String} containing the codon change string for the given alleleic information.
     */
    private static String getCodonChangeStringForInsertionInFrameWithInFrameStartPosition(final String alignedRefAllele,
                                                                                          final String alignedAltAllele,
                                                                                          final String nextRefCodon,
                                                                                          final int alignedCodingSequenceAlleleStart,
                                                                                          final int alignedReferenceAlleleStop ) {
        final StringBuilder ref = new StringBuilder();
        final StringBuilder alt = new StringBuilder();

        // For insertions we need to grab the codon just following the end of the insertion reference and place
        // it on the end of our reference before we render it.  We must also ignore the first codon because of
        // the leading extra base in the VCF file.
        // TODO: This is really bad - we are tying our output to a specific input format.  FIX IT.

        // Create a new reference allele by removing the first codon
        // (because of the added leading base for insertions)
        // and appending the next codon to it:
        final String newAlignedRefAllele = alignedRefAllele.substring(3) + nextRefCodon;

        // Create a new alternate allele by removing the first codon
        // (because of the added leading base for insertions):
        final String newAlignedAltAllele = alignedAltAllele.substring(3);

        // Create a new start and stop position because of the added leading base:
        final int newStartPosition = alignedCodingSequenceAlleleStart + 3;
        final int newEndPosition = alignedReferenceAlleleStop + 3;

        // Capitalize for insertion:

        // Insertions alternate alleles have ref first, then the extra bases,
        // so we capitalize only the extra bases.

        // The next codon is really our new ref allele.
        // And the ref will always be all lower case:
        ref.append(nextRefCodon.toLowerCase());

        // Add alleles to our stringstreams:
        // Alt allele chars will be upper case:
        alt.append(newAlignedAltAllele.toUpperCase());

        // Ref allele will be lower case:
        alt.append(newAlignedRefAllele.toLowerCase());

        // Here we add 3 to the start and stop positions because of the
        return "c.(" + newStartPosition + "-" + newEndPosition + ")" +
                ref.toString() + ">" + alt.toString();
    }


    /**
     * Get the codon change string for an in-frame Insertion with an out-of-frame start position.
     * @param alignedRefAllele The {@link String} of bases contained in the Reference allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param alignedAltAllele The {@link String} of bases contained in the Alternate allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param altAlleleLength The raw length of the Alternate allele as found in the input file (does not count any extra bases added to make the Alternate allele in-frame).
     * @param alignedCodingSequenceAlleleStart The position (1-based, inclusive) of the start of the first of the alleles in the coding sequence (prepended with up to 2 bases to be in-frame).
     * @param alignedReferenceAlleleStop The position (1-based, inclusive) of the in-frame end of the Reference allele.
     * @return The {@link String} containing the codon change string for the given alleleic information.
     */
    private static String getCodonChangeStringForInsertionInFrameWithOutOfFrameStartPosition(final String alignedRefAllele,
                                                                                             final String alignedAltAllele,
                                                                                             final int altAlleleLength,
                                                                                             final int alignedCodingSequenceAlleleStart,
                                                                                             final int alignedReferenceAlleleStop) {

        final StringBuilder alt = new StringBuilder();

        // For positions that are not in frame, we can just keep the ref and alt alleles since those will
        // correctly reflect the change:

        // Capitalize the insertion correctly:
        int i = 0;
        int j = 0;
        while ( alt.length() != alignedAltAllele.length() ) {
            if ( i == alignedRefAllele.length() ) {
                alt.append( alignedRefAllele.substring(i).toLowerCase() );
            }
            else {
                if ( alignedRefAllele.charAt(i) != alignedAltAllele.charAt(j) ) {
                    // Once we get 1 difference, we append the rest of the characters here:

                    // Get the rest of the characters from the alternate allele into the buffer:
                    for ( int k = 0 ; k < altAlleleLength - 1 ; ++k){
                        alt.append(Character.toUpperCase(alignedAltAllele.charAt(j++)));
                    }
                }
                else {
                    alt.append( Character.toLowerCase(alignedAltAllele.charAt(j)) );
                    ++i;
                    ++j;
                }
            }
        }

        return "c.(" + alignedCodingSequenceAlleleStart + "-" + alignedReferenceAlleleStop + ")" +
                alignedRefAllele.toLowerCase() + ">" + alt.toString();
    }

    /**
     * Get the codon change string for a Deletion.
     * @param alignedRefAllele The {@link String} of bases contained in the Reference allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param alignedAltAllele The {@link String} of bases contained in the Alternate allele (prepended and appended with up to 2 bases from the reference to be in-frame).
     * @param alignedCodingSequenceAlleleStart The position (1-based, inclusive) of the start of the first of the alleles in the coding sequence (prepended with up to 2 bases to be in-frame).
     * @param alignedReferenceAlleleStop The position (1-based, inclusive) of the in-frame end of the Reference allele.
     * @return The {@link String} containing the codon change string for the given alleleic information.
     */
    private static String getCodonChangeStringForDeletion(final String alignedRefAllele,
                                                          final String alignedAltAllele,
                                                          final int alignedCodingSequenceAlleleStart,
                                                          final int alignedReferenceAlleleStop) {
        // Check for frame shift syntax and return our string.
        // NOTE: We know that because of the way deletions are created
        // (i.e. a base is added to the beginning of the string), we can ignore the first codon (3 bases) of the
        // reference:
        String decorator = "del";
        if ( GATKProtectedVariantContextUtils.isFrameshift(alignedRefAllele, alignedAltAllele) ) {
            decorator = "fs";
        }

        return "c.(" + (alignedCodingSequenceAlleleStart + 3) + "-" +
                alignedReferenceAlleleStop + ")" +
                alignedRefAllele.substring(3).toLowerCase() + decorator;
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
     * Gets a codon change string for a splice site.
     * See {@link #createSpliceSiteCodonChange(int, int, int, int, Strand, int)} for details.
     */
    public static String createSpliceSiteCodonChange(final int variantStart,
                                                     final int exonNumber,
                                                     final int exonStart,
                                                     final int exonEnd,
                                                     final Strand strand) {
        return createSpliceSiteCodonChange(variantStart, exonNumber, exonStart, exonEnd, strand, 0);
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
     *     referenceCodingSequence
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
        if ( GATKProtectedVariantContextUtils.isOnp(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {
            return getProteinChangeStringForOnp(refAaSeq, altAaSeq, protChangeStartPos, protChangeEndPos);
        }
        // Check for the Insertion case:
        else if ( GATKProtectedVariantContextUtils.isInsertion(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {

            Utils.nonNull(seqComp.getCodingSequenceAlleleStart());
            Utils.nonNull(seqComp.getReferenceCodingSequence());
            assertValidStrand(seqComp.getStrand());

            // Because we have to deal with the codon after the insertion, we need to know what that next codon is:
            final String nextRefCodon = getNextReferenceCodon(seqComp.getReferenceCodingSequence(), seqComp.getAlignedCodingSequenceAlleleStart(), seqComp.getAlignedReferenceAlleleStop(), seqComp.getStrand());

            // We must also know the Amino Acid sequence for it:
            final String nextRefAaSeq = createAminoAcidSequence(nextRefCodon);

            if ( GATKProtectedVariantContextUtils.isFrameshift(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {
                return getProteinChangeStringForInsertionFrameshift(
                        seqComp.getCodingSequenceAlleleStart(),
                        refAaSeq,
                        protChangeStartPos,
                        nextRefAaSeq
                );
            }
            else {

                // If the variant position is in frame, then we have to use the next amino acid, because it will be the first
                // amino acid that is affected.
                // We add 1 because of the convention that the variant occurs just after the base specified for insertions.
                // TODO: This is really bad - we are tying our output to a specific input format.  FIX IT.
                if ( isPositionInFrame(seqComp.getCodingSequenceAlleleStart() + 1) ) {
                    return getProteinChangeStringForInsertionInFrameWithInFrameStartPosition(protChangeStartPos, altAaSeq);
                }
                else {
                    return getProteinChangeStringForInsertionInFrameWithOutOfFrameStartPosition(
                        refAaSeq,
                        altAaSeq,
                        protChangeStartPos,
                        seqComp.getAlternateAllele().length()
                    );
                }
            }
        }
        // Must be a deletion:
        else {
            return getProteinChangeStringForDeletion(seqComp.getReferenceAllele(),
                    seqComp.getAlternateAllele(),
                    seqComp.getReferenceAminoAcidSequence(),
                    seqComp.getProteinChangeStartPosition());
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
        if ( GATKProtectedVariantContextUtils.isFrameshift(referenceAllele, alternateAllele) ) {
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
        if ( GATKProtectedVariantContextUtils.isOnp(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {
            if (seqComp.getAlternateAllele().length() > 1) {
                return "c." + seqComp.getCodingSequenceAlleleStart() + "_" + (seqComp.getCodingSequenceAlleleStart() + seqComp.getReferenceAllele().length() - 1) +
                        seqComp.getReferenceAllele() + ">" + seqComp.getAlternateAllele();
            } else {
                return "c." + seqComp.getCodingSequenceAlleleStart() +
                        seqComp.getReferenceAllele() + ">" + seqComp.getAlternateAllele();
            }
        }
        // Check for Insertion:
        else if ( GATKProtectedVariantContextUtils.isInsertion(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {
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

            return "c." + start + "_" + end + "del" +
                    seqComp.getReferenceAllele().substring(seqComp.getAlternateAllele().length()).toUpperCase();
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
            alignedAlleleSeq = codingSequence.substring(start, end);
        }
        else {
            // Negative strand means we need to reverse complement and go from the other end:
            start = codingSequence.length() - alignedAlleleStop;
            end = codingSequence.length() - alignedAlleleStart + 1;
            alignedAlleleSeq = ReadUtils.getBasesReverseComplement( codingSequence.substring(start, end).getBytes() );
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

        // Set the window on our reference to be correct for our start and end:
        reference.setWindow(
                Math.abs(start - reference.getInterval().getStart()),
                Math.abs(reference.getInterval().getEnd() - end)
        );

        // Now that the window size is correct, we can go through and pull our sequences out.

        // Get the window so we can convert to reference coordinates from genomic coordinates of the exons:
        final SimpleInterval refWindow = reference.getWindow();
        final byte[] bases = reference.getBases();

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

    /**
     * A simple data object to hold a comparison between a reference sequence and an alternate allele.
     */
    public static class SequenceComparison {

        /**
         * Bases covering the region around the variant.
         * This includes bases that would overlap any of the variant or reference bases.
         * (i.e. if the reference allele is 'A' and the variant is 'GCGCG', this would be
         * 'ANNNNN' where 'N' is the correct base at that position in the reference sequence).
         * Contains {@link FuncotatorUtils.SequenceComparison#referenceWindow} bases of padding before and
         * after the total length.
         */
        private String referenceBases                     = null;

        /**
         * The number of bases in {@link FuncotatorUtils.SequenceComparison#referenceBases} before the
         * start of the reference Allele / variant.
         */
        private Integer referenceWindow                   = null;

        /**
         * The reference coding sequence for a the transcript of this sequence comparison.
         * This does NOT include introns.
         * Stored in the forward reading direction.  For NEGATIVE strand reads, must
         * reverse complement any bases retrieved.
         */
        private ReferenceSequence referenceCodingSequence = null;

        /**
         * The contig on which this sequence comparison occurs.
         */
        private String  contig                               = null;

        /**
         * The strand on which this sequence comparison occurs.
         */
        private Strand  strand                               = null;

        /**
         * The position (1-based, inclusive in genome coordinates - relative to the start of
         * {@link FuncotatorUtils.SequenceComparison#contig}) of the start of the reference allele / variant.
         */
        private Integer alleleStart                          = null;

        /**
         * The position (1-based, inclusive) in transcript coordinates relative to the start of
         * the transcript of start of the reference allele / variant.
         */
        private Integer transcriptAlleleStart                = null;

        /**
         * The position (1-based, inclusive) in coding sequence coordinates relative to the start of
         * the coding region of the transcript of the start of the reference allele / variant.
         * This location is obtained by concatenating the exons together and counting from the start of
         * that sequence to where the reference allele / variant begins.
         */
        private Integer codingSequenceAlleleStart            = null;

        /**
         * The in-frame position (1-based, inclusive) in coding sequence coordinates relative to the start of
         * the coding region of the transcript of the start of the first codon containing the reference allele / variant.
         * This location is obtained by concatenating the exons together and counting from the start of
         * that sequence to where the reference allele / variant begins, then moving backwards to an in-frame
         * position, if necessary.
         */
        private Integer alignedCodingSequenceAlleleStart     = null;

        /**
         * The position (1-based, inclusive in genome coordinates - relative to the start of
         * {@link FuncotatorUtils.SequenceComparison#contig}) of the start of the exon that contains the
         * variant in this {@link FuncotatorUtils.SequenceComparison}.
         */
        private Integer exonStartPosition                    = null;

        /**
         * The position (1-based, inclusive in genome coordinates - relative to the start of
         * {@link FuncotatorUtils.SequenceComparison#contig}) of the end of the exon that contains the
         * variant in this {@link FuncotatorUtils.SequenceComparison}.
         */
        private Integer exonEndPosition                      = null;

        /**
         * The start position (1-based, inclusive) of the Protein Change for the alleles in this {@link SequenceComparison}.
         * This is computed using {@link SequenceComparison#alignedCodingSequenceAlleleStart}.
         * See {@link FuncotatorUtils#getProteinChangePosition} for more details.
         */
        private Integer proteinChangeStartPosition           = null;

        /**
         * The end position (1-based, inclusive) of the Protein Change for the alleles in this {@link SequenceComparison}.
         * This is computed by using {@link SequenceComparison#alignedCodingSequenceAlleleStart} and the length of
         * {@link SequenceComparison#alignedCodingSequenceAlternateAllele}.
         * See {@link FuncotatorUtils#getProteinChangeEndPosition} for more details.
         */
        private Integer proteinChangeEndPosition             = null;

        /**
         * A string representation of the reference allele.
         */
        private String  referenceAllele                      = null;
        /**
         * An in-frame sequence of bases that overlaps the given reference allele based on the raw reference genome.
         * (i.e. This includes INTRONS.)
         */
        private String  alignedReferenceAllele               = null;
        /**
         * An in-frame sequence of bases that overlaps the given reference allele for the coding region only.
         * (i.e. This includes ONLY EXONS.)
         */
        private String  alignedCodingSequenceReferenceAllele = null;

        /**
         * The in-frame position (1-based, inclusive) of the last base of the last codon that contains the reference
         * allele relative to the start of the coding sequence.
         * All codons containing the reference allele can be extracted from the coding sequence using
         * {@link SequenceComparison#alignedCodingSequenceAlleleStart} and {@link SequenceComparison#alignedReferenceAlleleStop}.
         */
        private Integer alignedReferenceAlleleStop           = null;

        /**
         * The amino acid sequence as coded by the in-frame reference coding sequence ({@link SequenceComparison#alignedCodingSequenceReferenceAllele}.
         */
        private String  referenceAminoAcidSequence           = null;

        /**
         * A string representation of the alternate allele.
         */
        private String  alternateAllele                      = null;

        /**
         * An in-frame sequence of bases that includes the entire alternate allele based on the reference genome.
         * May span multiple codons.
         * May include intron bases.
         */
        private String  alignedAlternateAllele               = null;

        /**
         * An in-frame sequence of bases that includes the entire alternate allele based on the coding sequence reference.
         * May span multiple codons.
         */
        private String  alignedCodingSequenceAlternateAllele = null;

        /**
         * The in-frame position (1-based, inclusive) of the last base of the last codon that contains the alternate
         * allele relative to the start of the coding sequence.
         * All codons containing the alternate allele can be extracted from the coding sequence using
         * {@link SequenceComparison#alignedCodingSequenceAlleleStart} and {@link SequenceComparison#alignedAlternateAlleleStop}.
         */
        private Integer alignedAlternateAlleleStop           = null;

        /**
         * The amino acid sequence as coded by the in-frame alternate coding sequence ({@link SequenceComparison#alignedCodingSequenceAlternateAllele}.
         */
        private String  alternateAminoAcidSequence           = null;

        // =============================================================================================================

        public String getReferenceBases() {
            return referenceBases;
        }

        public void setReferenceBases(final String referenceBases) {
            this.referenceBases = referenceBases;
        }

        public Integer getReferenceWindow() {
            return referenceWindow;
        }

        public void setReferenceWindow(final Integer referenceWindow) {
            this.referenceWindow = referenceWindow;
        }

        /**
         * Return the {@link ReferenceSequence} containing the coding region for the transcript of this {@link SequenceComparison}.
         * This does NOT include introns.
         * The reference sequence is stored in the forward reading direction.
         * For NEGATIVE strand reads, must reverse complement any bases retrieved.
         */
        public ReferenceSequence getReferenceCodingSequence() {
            return referenceCodingSequence;
        }

        public void setReferenceCodingSequence(final ReferenceSequence referenceCodingSequence) {
            this.referenceCodingSequence = referenceCodingSequence;
        }

        public String getContig() {
            return contig;
        }

        public void setContig(final String contig) {
            this.contig = contig;
        }

        public Strand getStrand() {
            return strand;
        }

        public void setStrand(final Strand strand) {

            if (strand == Strand.NONE) {
                throw new GATKException("Cannot handle NONE strand.");
            }

            this.strand = strand;
        }

        public Integer getAlleleStart() {
            return alleleStart;
        }

        public void setAlleleStart(final Integer alleleStart) {
            this.alleleStart = alleleStart;
        }

        public Integer getTranscriptAlleleStart() {
            return transcriptAlleleStart;
        }

        public void setTranscriptAlleleStart(final Integer transcriptAlleleStart) {
            this.transcriptAlleleStart = transcriptAlleleStart;
        }

        public Integer getCodingSequenceAlleleStart() {
            return codingSequenceAlleleStart;
        }

        public void setCodingSequenceAlleleStart(final Integer codingSequenceAlleleStart) {
            this.codingSequenceAlleleStart = codingSequenceAlleleStart;
        }

        public Integer getAlignedCodingSequenceAlleleStart() {
            return alignedCodingSequenceAlleleStart;
        }

        public void setAlignedCodingSequenceAlleleStart(final Integer alignedCodingSequenceAlleleStart) {
            this.alignedCodingSequenceAlleleStart = alignedCodingSequenceAlleleStart;
        }

        public Integer getExonStartPosition() {
            return exonStartPosition;
        }

        public void setExonStartPosition(final Integer exonStartPosition) {
            this.exonStartPosition = exonStartPosition;
        }

        public Integer getExonEndPosition() {
            return exonEndPosition;
        }

        public void setExonEndPosition(final Integer exonEndPosition) {
            this.exonEndPosition = exonEndPosition;
        }

        public void setExonPosition( final SimpleInterval exonPosition ) {
            this.exonStartPosition = exonPosition.getStart();
            this.exonEndPosition = exonPosition.getEnd();
        }

        public Integer getProteinChangeStartPosition() {
            return proteinChangeStartPosition;
        }

        public void setProteinChangeStartPosition(final Integer proteinChangeStartPosition) {
            this.proteinChangeStartPosition = proteinChangeStartPosition;
        }

        public Integer getProteinChangeEndPosition() {
            return proteinChangeEndPosition;
        }

        public void setProteinChangeEndPosition(final Integer proteinChangeEndPosition) {
            this.proteinChangeEndPosition = proteinChangeEndPosition;
        }

        public String getReferenceAllele() {
            return referenceAllele;
        }

        public void setReferenceAllele(final String referenceAllele) {
            this.referenceAllele = referenceAllele;
        }

        public String getAlignedReferenceAllele() {
            return alignedReferenceAllele;
        }

        public void setAlignedReferenceAllele(final String alignedReferenceAllele) {
            this.alignedReferenceAllele = alignedReferenceAllele;
        }

        public String getAlignedCodingSequenceReferenceAllele() {
            return alignedCodingSequenceReferenceAllele;
        }

        public void setAlignedCodingSequenceReferenceAllele(final String alignedCodingSequenceReferenceAllele) {
            this.alignedCodingSequenceReferenceAllele = alignedCodingSequenceReferenceAllele;
        }

        public Integer getAlignedReferenceAlleleStop() {
            return alignedReferenceAlleleStop;
        }

        public void setAlignedReferenceAlleleStop(final Integer alignedReferenceAlleleStop) {
            this.alignedReferenceAlleleStop = alignedReferenceAlleleStop;
        }

        public String getReferenceAminoAcidSequence() {
            return referenceAminoAcidSequence;
        }

        public void setReferenceAminoAcidSequence(final String referenceAminoAcidSequence) {
            this.referenceAminoAcidSequence = referenceAminoAcidSequence;
        }

        public String getAlternateAllele() {
            return alternateAllele;
        }

        public void setAlternateAllele(final String alternateAllele) {
            this.alternateAllele = alternateAllele;
        }

        public String getAlignedAlternateAllele() {
            return alignedAlternateAllele;
        }

        public void setAlignedAlternateAllele(final String alignedAlternateAllele) {
            this.alignedAlternateAllele = alignedAlternateAllele;
        }

        public String getAlignedCodingSequenceAlternateAllele() {
            return alignedCodingSequenceAlternateAllele;
        }

        public void setAlignedCodingSequenceAlternateAllele(final String alignedCodingSequenceAlternateAllele) {
            this.alignedCodingSequenceAlternateAllele = alignedCodingSequenceAlternateAllele;
        }

        public Integer getAlignedAlternateAlleleStop() {
            return alignedAlternateAlleleStop;
        }

        public void setAlignedAlternateAlleleStop(final Integer alignedAlternateAlleleStop) {
            this.alignedAlternateAlleleStop = alignedAlternateAlleleStop;
        }

        public String getAlternateAminoAcidSequence() {
            return alternateAminoAcidSequence;
        }

        public void setAlternateAminoAcidSequence(final String alternateAminoAcidSequence) {
            this.alternateAminoAcidSequence = alternateAminoAcidSequence;
        }
    }
}
