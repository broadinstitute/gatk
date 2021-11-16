package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.builder.fluent.Configurations;
import org.apache.commons.configuration2.ex.ConfigurationException;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public final class FuncotatorUtils {

    private final static Logger logger = LogManager.getLogger(FuncotatorUtils.class);

    /**
     * PRIVATE CONSTRUCTOR
     * DO NOT INSTANTIATE THIS CLASS!
     */
    private FuncotatorUtils() {}

    public static final int DEFAULT_MIN_NUM_BASES_FOR_VALID_SEGMENT = 150;

    private static final Map<String, AminoAcid> tableByCodon;
    private static final Map<String, AminoAcid> tableByCode;
    private static final Map<String, AminoAcid> tableByLetter;

    private static final Map<String, AminoAcid>              mtDifferentAaTableByCodon;
    private static final Map<Pair<Genus, String>, AminoAcid> mtSpecialStartCodonsBySpecies;

    private static SAMSequenceDictionary B37_SEQUENCE_DICTIONARY = null;

    private static final Map<String, String> B37_To_HG19_CONTIG_NAME_MAP;
    private static final Map<String, String> HG19_TO_B37_CONTIG_NAME_MAP;

    private static final String CODON_CHANGE_FORMAT_STRING   = "c.(%d-%d)%s";
    private static final String PROTEIN_CHANGE_FORMAT_STRING = "p.%s%s%s%s";
    private static final String CDNA_CHANGE_FORMAT_STRING    = "c.%s%s%s%s";

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

        final HashMap<String, AminoAcid>              mtCodons        = new HashMap<>();
        final HashMap<Pair<Genus, String>, AminoAcid> mtSpecialStarts = new HashMap<>();

        // Add in the codons that are different in the MT code, including IUPAC bases:
        mtCodons.put("ATA", AminoAcid.METHIONINE);
        mtCodons.put("AGA", AminoAcid.STOP_CODON);
        mtCodons.put("AGG", AminoAcid.STOP_CODON);
        mtCodons.put("AGR", AminoAcid.STOP_CODON);
        mtCodons.put("TGA", AminoAcid.TRYPTOPHAN);

        // Add in the codons that serve as alternate start codons in various organisms, including IUPAC bases:
        mtSpecialStarts.put(Pair.of(Genus.BOS, "ATA"), AminoAcid.METHIONINE);
        mtSpecialStarts.put(Pair.of(Genus.HOMO, "ATT"), AminoAcid.METHIONINE);
        mtSpecialStarts.put(Pair.of(Genus.MUS, "ATT"), AminoAcid.METHIONINE);
        mtSpecialStarts.put(Pair.of(Genus.MUS, "ATC"), AminoAcid.METHIONINE);
        mtSpecialStarts.put(Pair.of(Genus.MUS, "ATY"), AminoAcid.METHIONINE);
        mtSpecialStarts.put(Pair.of(Genus.CORTURNIX, "GTG"), AminoAcid.METHIONINE);
        mtSpecialStarts.put(Pair.of(Genus.GALLUS, "GTG"), AminoAcid.METHIONINE);


        tableByCodon = Collections.unmodifiableMap(mapByCodon);
        tableByCode = Collections.unmodifiableMap(mapByCode);
        tableByLetter = Collections.unmodifiableMap(mapByLetter);
        mtDifferentAaTableByCodon = Collections.unmodifiableMap(mtCodons);
        mtSpecialStartCodonsBySpecies = Collections.unmodifiableMap(mtSpecialStarts);

        B37_To_HG19_CONTIG_NAME_MAP = initializeB37ToHg19ContigNameMap();
        HG19_TO_B37_CONTIG_NAME_MAP = B37_To_HG19_CONTIG_NAME_MAP.entrySet()
                                            .stream()
                                            .collect(Collectors.toMap(Map.Entry::getValue, Map.Entry::getKey));
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
     * Returns the {@link AminoAcid} corresponding to the given three-letter Mitochondrial {@code rawCodon}.
     * The codons given are expected to be valid for Mitochondrial DNA.
     * Assumes no special cases for alternate initiation sites.
     * @param rawCodon The three-letter codon (each letter one of A,[T or U],G,C) representing a Mitochondrial {@link AminoAcid}.  Must not be {@code null}.
     * @return The {@link AminoAcid} corresponding to the given {@code rawCodon}.  Returns {@code null} if the given {@code rawCodon} does not code for a Mitochondrial {@link AminoAcid}.
     */
    public static AminoAcid getMitochondrialAminoAcidByCodon(final String rawCodon) {
        return getMitochondrialAminoAcidByCodon(rawCodon, false, Genus.UNSPECIFIED);
    }

    /**
     * Returns the {@link AminoAcid} corresponding to the given three-letter Mitochondrial {@code rawCodon}.
     * The codons given are expected to be valid for Mitochondrial DNA.
     * @param rawCodon The three-letter codon (each letter one of A,[T or U],G,C) representing a Mitochondrial {@link AminoAcid}
     * @param isFirst {@code true} iff the given codon appears first in the coding sequence.  {@code false} otherwise.
     * @return The {@link AminoAcid} corresponding to the given {@code rawCodon}.  Returns {@code null} if the given {@code rawCodon} does not code for a Mitochondrial {@link AminoAcid}.
     */
    public static AminoAcid getMitochondrialAminoAcidByCodon(final String rawCodon,
                                                             final boolean isFirst,
                                                             final Genus genus) {
        if (rawCodon == null) {
            return null;
        }

        // TODO: Need to solicit more info on partly coded stop codons as alluded to here: https://www.sciencedirect.com/science/article/pii/S0005272898001613 (Issue https://github.com/broadinstitute/gatk/issues/5363)

        // Convert Uracils to Thymines and convert to upper case so we can use our normal lookup table:
        // Note: this may be unnecessary, but is here for safety and correctness (at least according to the
        //       lookup table on wikipedia https://en.wikipedia.org/wiki/Vertebrate_mitochondrial_code).
        final String upperCodon = rawCodon.replaceAll("[Uu]", "T").toUpperCase();

        // Check special species-specific start codons here:
        if ( isFirst && mtSpecialStartCodonsBySpecies.containsKey(Pair.of(genus, upperCodon)) ) {
            return mtSpecialStartCodonsBySpecies.get(Pair.of(genus, upperCodon));
        }
        // Check for MT contig-specific codons here:
        else if ( mtDifferentAaTableByCodon.containsKey(upperCodon) ) {
            return mtDifferentAaTableByCodon.get(upperCodon);
        }
        // Everything else is the same as the Standard Code:
        else {
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

        return (((position - 1) % AminoAcid.CODON_LENGTH) == 0);
    }

    /**
     * Get the start position of the difference between the reference and alternate alleles in the given {@link VariantContext}.
     * Assumes the given {@link VariantContext} contains only one alternate allele.
     * @param variant {@link VariantContext} in which to determine the start of the different bases.
     * @return The start position of the difference between the reference and alternate alleles in the given {@link VariantContext}.
     */
    public static int getIndelAdjustedAlleleChangeStartPosition(final VariantContext variant) {
        return getIndelAdjustedAlleleChangeStartPosition(variant, variant.getAlternateAllele(0));
    }

    /**
     * Get the start position of the difference between the reference and alternate alleles in the given {@link VariantContext}.
     * @param variant {@link VariantContext} in which to determine the start of the different bases.
     * @param altAllele The alternate {@link Allele} against which to check the reference allele in {@code variant}.
     * @return The start position of the difference between the reference and alternate alleles in the given {@link VariantContext}.
     */
    public static int getIndelAdjustedAlleleChangeStartPosition(final VariantContext variant, final Allele altAllele) {

        // If the variant is an indel, we need to check only the bases that are added/deleted for overlap.
        // The convention for alleles in Funcotator is to preserve a leading base for an indel, so we just need
        // to create a new variant that has its start position shifted by the leading base.
        // NOTE: because there could be degenerate VCF files that have more than one leading base overlapping, we need
        //       to detect how many leading bases there are that overlap, rather than assuming there is only one.
        final int varStart;
        if ( GATKVariantContextUtils.typeOfVariant(variant.getReference(), altAllele).equals(VariantContext.Type.INDEL) &&
             !GATKVariantContextUtils.isComplexIndel(variant.getReference(), altAllele) ) {
            int startOffset = 0;
            while ( (startOffset < variant.getReference().length()) && (startOffset < altAllele.length()) && (variant.getReference().getBases()[ startOffset ] == altAllele.getBases()[ startOffset ]) ) {
                ++startOffset;
            }
            varStart = variant.getStart() + startOffset;
        }
        else {
            // Not an indel?  Then we should have no overlapping bases:
            varStart = variant.getStart();
        }

        return varStart;
    }

    /**
     * Get the string of bases that are different from the given alleles.
     * Assumes there is one contiguous string of changed bases between the two alleles.
     * Assumes that if there is overlap between the alleles, the overlap occurs at either the front or the back.
     * @param firstAllele First {@link Allele}.  Must not be {@code null}.
     * @param secondAllele Second {@link Allele}.  Must not be {@code null}.
     * @param copyRefBasesWhenAltIsPastEnd Will copy the bases from the given {@code firstAllele} when the alternate allele has no more bases to copy over.  Used primarily for handling deletions.
     * @return A string containing the bases from the given {@code secondAllele} that are different from the {@code firstAllele} (in their correct relative order).
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
    public static int getStartPositionInTranscript(final Locatable variant,
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

        return (int)(Math.ceil(alleleEndPosition / ((double)AminoAcid.CODON_LENGTH)) * AminoAcid.CODON_LENGTH);
    }

    /**
     * Gets the sequence aligned position (1-based, inclusive) for the given coding sequence position.
     * This will produce the next lowest position evenly divisible by {@link AminoAcid#CODON_LENGTH}, such that a codon starting at this returned
     * position would include the given position.  This can be a negative number, in which case the codon would start
     * at the given position relative to the normal starting position (1) (this would be the case for an upstream
     * UTR or flank).
     *
     * @param position A sequence starting coordinate for which to produce an coding-aligned position.
     * @return A coding-aligned position (1-based, inclusive) corresponding to the given {@code position}.
     */
    public static int getAlignedPosition(final int position) {

        if ( position > 0 ) {
            return position - ((position - 1) % AminoAcid.CODON_LENGTH);
        }
        else {
            final int adjustedPos = 1 - position;
            return -(adjustedPos - ((adjustedPos - 1) % AminoAcid.CODON_LENGTH) + 1);
        }
    }

    /**
     * Calculates whether the given {@code startPosition} (1-based, inclusive) is in frame relative to the end of the region.
     * @param startPosition The position (1-based, inclusive) relative to the start of a region to check for frame alignment.  Must be relative to the given {@code regionLength}.  Concretely, {@code 0} < {@code startPosition} <= {@code regionLength}.
     * @param regionLength The length of the region containing {@code startPosition}.  Must be >= 0.
     * @return {@code true} if the given {@code startPosition} is in frame relative to the given {@code regionLength} ; {@code false} otherwise.
     */
    public static boolean isInFrameWithEndOfRegion(final int startPosition, final int regionLength) {

        ParamUtils.isPositiveOrZero( regionLength, "Region length must be >= 0." );

        // Micro-optimization to split the return statements based on
        // the if statement so we don't have to do unnecessary math in the "normal" case:
        if ( startPosition > 0 ) {
            // Add 1 because of the 1-based/inclusive positions:
            return (((regionLength - startPosition + 1) % AminoAcid.CODON_LENGTH) == 0);
        }
        else {
            // If we have a position before our region starts, we must calculate an offset
            // between the position and the actual start.
            // We can then add this offset to the start and region length to simplify the calculation.
            final int preFlankOffset = 1 - startPosition;

            return ((((regionLength + preFlankOffset) - (startPosition + preFlankOffset) + 1) % AminoAcid.CODON_LENGTH) == 0);
        }
    }

    /**
     * Checks to see whether a given indel location occurs on a codon boundary.
     * That is, whether the given indel location is not within a codon but is cleanly between two adjacent codons.
     * NOTE: ASSUMES that there is a leading base that is not part of the indel prepended to the indel string for context.
     * @param codingSequenceAlleleStart The start position of the variant in the coding sequence.
     * @param alignedCodingSequenceAlleleStart The start position of the first codon containing part of the variant in the coding sequence.
     * @param refAllele A {@link String} containing the bases in the reference allele for the given variant. Must not be {@code null}.
     * @param strand The {@link Strand} on which the current codon resides.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return {@code true} if the given indel cleanly occurs between two adjacent codons; {@code false} otherwise.
     */
    public static boolean isIndelBetweenCodons(final int codingSequenceAlleleStart,
                                                final int alignedCodingSequenceAlleleStart,
                                                final String refAllele,
                                                final Strand strand) {
        Utils.nonNull(refAllele);
        assertValidStrand(strand);

        if ( strand == Strand.POSITIVE ) {
            // Check normally for positive strands:
            final int codonOffset = codingSequenceAlleleStart - alignedCodingSequenceAlleleStart;
            return (((codonOffset + refAllele.length()) % AminoAcid.CODON_LENGTH) == 0);
        }
        else {
            // This is a little strange. Because we're on the reverse strand the bases will be inserted in
            // the opposite direction than what we expect (even though the coordinates are now in transcription order).
            // Because of this, we just need to check if we're at the start of a codon for if this occurs between them.
            //
            // For example:
            //
            //   1:1249193 C->CGCAA
            //
            //           insertion here
            //                 |
            //                 V
            //   +   ...aa|agt|ctt|gcg|ga...
            //   -   ...tt|tca|gaa|cgc|ct...
            //
            //                 |
            //                  \___
            //                  |   |
            //                  V   V
            //   +   ...aaa|gtc|GCA|Att|gcg|ga...
            //   -   ...ttt|cag|CGT|Taa|cgc|ct...
            //
            // If we inserted them on the + strand the bases would be inserted between codons, but because we're on
            // the - strand, the bases are inserted within the `aag` codon.

            return (codingSequenceAlleleStart == alignedCodingSequenceAlleleStart);
        }
    }

    /**
     * Create a properly capitalized version of the given alternate allele codon change string for an insertion by
     * capitalizing only the bases that are different from the reference.
     * This is done using the index of where the alternate allele starts and the length of the alternate allele.
     * @param alignedAlternateCodingSequenceAllele The alternate allele codon {@link String} to capitalize properly.
     * @param alternateAllele A {@link String} containing the bases in the alternate allele.
     * @param startingOffset The offset in {@code alternateAlleleCodonChangeString} where the alternate allele begins.
     * @return A properly capitalized version of the given {@code alternateAlleleCodonChangeString} based on the given alternate allele information.
     */
    private static String createCapitalizedAlternateAlleleInsertionCodonChangeString(final String alignedAlternateCodingSequenceAllele,
                                                                                     final String alternateAllele,
                                                                                     final int startingOffset ) {
        final StringBuilder sb = new StringBuilder();

        // Account for the leading base that we require for insertions:
        final int newStartingOffset = startingOffset + 1;

        int i = 0;
        while ( i < newStartingOffset ) {
            sb.append(Character.toLowerCase( alignedAlternateCodingSequenceAllele.charAt(i++) ) );
        }
        // Subtract 1 from the allele length to account for the leading base that we require for insertions:
        while ( (i - newStartingOffset) < (alternateAllele.length() - 1) ) {
            sb.append(Character.toUpperCase( alignedAlternateCodingSequenceAllele.charAt(i++) ) );
        }
        while ( i < alignedAlternateCodingSequenceAllele.length() ) {
            sb.append(Character.toLowerCase( alignedAlternateCodingSequenceAllele.charAt(i++) ) );
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
     *     strand
     * For INDELS, these are also required:
     *     contig,
     *     alleleStart
     *     transcriptCodingSequence
     *     alignedAlternateAllele
     * @param seqComp {@link SequenceComparison} representing the alternate and reference alleles for a DNA sequence.  Must not be {@code null}.
     * @param startCodon {@link Locatable} representing the start codon for the coding region for the variant in {@code seqComp}.  May be {@code null}.
     * @return A {@link String} representing the codon change for the given {@link SequenceComparison}.
     */
    public static String getCodonChangeString( final SequenceComparison seqComp,
                                               final Locatable startCodon ) {
        Utils.nonNull(seqComp.getAlignedCodingSequenceAlleleStart());
        Utils.nonNull(seqComp.getAlignedReferenceAlleleStop());
        Utils.nonNull(seqComp.getReferenceAllele());
        Utils.nonNull(seqComp.getAlignedCodingSequenceReferenceAllele());
        Utils.nonNull(seqComp.getAlternateAllele());
        Utils.nonNull(seqComp.getAlignedCodingSequenceAlternateAllele());
        Utils.nonNull(seqComp.getCodingSequenceAlleleStart());
        assertValidStrand(seqComp.getStrand());

        // WHY would seqComp.getAlignedCodingSequenceReferenceAllele() == seqComp.getAlignedCodingSequenceAlternateAllele() here?
        // see 1:1152971 T>C for details

        // ONP:
        if ( GATKVariantContextUtils.isXnp(seqComp.getAlignedCodingSequenceReferenceAllele(), seqComp.getAlignedCodingSequenceAlternateAllele()) ) {
            return getCodonChangeStringForOnp(
                    seqComp.getAlignedCodingSequenceReferenceAllele(),
                    seqComp.getAlignedCodingSequenceAlternateAllele(),
                    seqComp.getAlignedCodingSequenceAlleleStart(),
                    seqComp.getAlignedReferenceAlleleStop()
            );
        }
        // INDELS:
        else {
            // Get whether the indel occurs cleanly between codons:
            final boolean indelIsBetweenCodons =
                    isIndelBetweenCodons(seqComp.getCodingSequenceAlleleStart(),
                            seqComp.getAlignedCodingSequenceAlleleStart(),
                            seqComp.getReferenceAllele(),
                            seqComp.getStrand());

            // Cache whether this is an insertion or deletion.
            // Note: If the variant is not an insertion, it must be a deletion.
            final boolean isInsertion = GATKVariantContextUtils.isInsertion(seqComp.getAlignedCodingSequenceReferenceAllele(), seqComp.getAlignedCodingSequenceAlternateAllele());

            // Get our baseic aligned start and stop positions:
            final int alignedCodonStart = seqComp.getAlignedCodingSequenceAlleleStart();
            // Subtract 1 for the inclusive variant positions:
            final int alignedCodonEnd = seqComp.getAlignedCodingSequenceAlleleStart() + seqComp.getAlignedCodingSequenceReferenceAllele().length() - 1;

            // Get the ref codon for the positions where the indel is:
            final String refCodon = seqComp.getAlignedCodingSequenceReferenceAllele().toLowerCase();

            // If we have a start codon indel, we do not need to render the codon change string.
            // (This is an oncotator convention that has been carried over.)
            if ( isStartCodonIndel(seqComp, startCodon) ) {
                return "";
            }
            // Frameshifts are all the same:
            else if ( GATKVariantContextUtils.isFrameshift(seqComp.getAlignedCodingSequenceReferenceAllele(), seqComp.getAlignedCodingSequenceAlternateAllele()) ) {
                return getCodonChangeStringForFrameShifts(seqComp, isInsertion, alignedCodonStart, alignedCodonEnd, refCodon);
            }
            else {
                // Insertions:
                if ( isInsertion ) {
                    return getCodonChangeStringForInsertion(seqComp, indelIsBetweenCodons, alignedCodonStart, alignedCodonEnd, refCodon);
                }
                // Deletions:
                else {
                    return getCodonChangeStringForDeletion(seqComp, indelIsBetweenCodons, alignedCodonStart, alignedCodonEnd, refCodon);
                }
            }
        }
    }

    private static String getCodonChangeStringForDeletion(final SequenceComparison seqComp, final boolean indelIsBetweenCodons, int alignedCodonStart, int alignedCodonEnd, String refCodon) {

        // Requires:
        //     seqComp.getStrand()
        //     seqComp.getAlignedCodingSequenceAlternateAllele()

        if ( indelIsBetweenCodons ) {

            if ( seqComp.getStrand() == Strand.POSITIVE ) {
                // Skip the first AminoAcid.CODON_LENGTH bases in the aligned codon:
                alignedCodonStart += AminoAcid.CODON_LENGTH;
                refCodon = refCodon.substring(AminoAcid.CODON_LENGTH);
            }
            else {
                // Skip the last AminoAcid.CODON_LENGTH bases in the aligned codon:
                alignedCodonEnd -= AminoAcid.CODON_LENGTH;
                refCodon = refCodon.substring(0, refCodon.length() - AminoAcid.CODON_LENGTH);
            }

            return String.format(CODON_CHANGE_FORMAT_STRING, alignedCodonStart, alignedCodonEnd, refCodon + "del");
        }
        else {
            return String.format(CODON_CHANGE_FORMAT_STRING, alignedCodonStart, alignedCodonEnd, refCodon + ">" + seqComp.getAlignedCodingSequenceAlternateAllele().toLowerCase());
        }
    }

    private static String getCodonChangeStringForInsertion(final SequenceComparison seqComp, final boolean indelIsBetweenCodons, int alignedCodonStart, final int alignedCodonEnd, final String refCodon) {

        // Requires:
        //     seqComp.getTranscriptCodingSequence()
        //     seqComp.getAlignedCodingSequenceAlleleStart()
        //     seqComp.getAlignedReferenceAlleleStop()
        //     seqComp.getStrand()
        //     seqComp.getAlignedCodingSequenceAlternateAllele()
        //     seqComp.getAlternateAllele()
        //     seqComp.getCodingSequenceAlleleStart()
        //     seqComp.getAlignedCodingSequenceReferenceAllele()
        //     seqComp.getAlignedAlternateAllele()
        //     seqComp.getStrand()
        //     seqComp.getStrand()

        if ( indelIsBetweenCodons ) {

            // Get the codon that is either after or before this position:
            final String adjacentRefCodon = getAdjacentReferenceCodon(
                    seqComp.getTranscriptCodingSequence(),
                    seqComp.getAlignedCodingSequenceAlleleStart(),
                    seqComp.getAlignedReferenceAlleleStop(),
                    seqComp.getStrand()
            );

            // If we're at the start or end of the sequence, we need to handle it as a special case:
            if ( adjacentRefCodon.isEmpty() ) {
                // We're at the very start/end of our coding sequence, so it doesn't make sense to include an
                // extra codon at the beginning/end (because there is no extra codon).

                // Get the capitalization correct in the alt string:
                final String capitalizedAltString = createCapitalizedAlternateAlleleInsertionCodonChangeString(
                        seqComp.getAlignedCodingSequenceAlternateAllele(),
                        seqComp.getAlternateAllele(),
                        seqComp.getCodingSequenceAlleleStart() - seqComp.getAlignedCodingSequenceAlleleStart() -
                                // Need to subtract 1 if we're going on the other strand to line up correctly with
                                // the start of the variant alleles:
                                (seqComp.getStrand() == Strand.POSITIVE ? 0 : 1)
                );

                return String.format(CODON_CHANGE_FORMAT_STRING,
                        seqComp.getAlignedCodingSequenceAlleleStart(),
                        (seqComp.getAlignedCodingSequenceAlleleStart() + seqComp.getAlignedCodingSequenceReferenceAllele().length()-1),
                        seqComp.getAlignedCodingSequenceReferenceAllele().toLowerCase() + ">" + capitalizedAltString);
            }
            // For + strand, we just get the changed codon, then the ref codon:
            else if (seqComp.getStrand() == Strand.POSITIVE ) {

                // in the + direction, the adjacent codon is the NEXT codon:
                final String nextRefCodon = adjacentRefCodon;

                // Here we adjust everything "right" by AminoAcid.CODON_LENGTH bases (1 codon) because of the leading base that is
                // required for indels in VCF:
                return String.format(CODON_CHANGE_FORMAT_STRING,
                        (seqComp.getAlignedCodingSequenceAlleleStart() + AminoAcid.CODON_LENGTH),
                        (seqComp.getAlignedCodingSequenceAlleleStart() + AminoAcid.CODON_LENGTH + 2),
                        nextRefCodon.toLowerCase() + ">" + seqComp.getAlignedAlternateAllele().substring(AminoAcid.CODON_LENGTH) + nextRefCodon.toLowerCase());
            }
            // For - strand we get both the codon before and after the insertion:
            else {
                // in the - direction, the adjacent codon is the PREVIOUS codon:
                final String prevRefCodon = adjacentRefCodon;

                // Adjust start for the extra codons:
                alignedCodonStart -= AminoAcid.CODON_LENGTH;

                // NOTE: The alt allele will contain the NEXT ref codon as well as the alt.
                //       The first codons in `seqComp.getAlignedAlternateAllele()` will be the alt, the last
                //       will be the ref.

                // Get the real alt codon:
                final String altCodon = seqComp.getAlignedAlternateAllele().substring(0, seqComp.getAlignedAlternateAllele().length() - AminoAcid.CODON_LENGTH);

                // NOTE: `nextRefCodon` will actually contain the previous ref codon because we're on the -
                //       strand.  Yes, it's confusing.  Sorry about that.
                return String.format(CODON_CHANGE_FORMAT_STRING,
                        alignedCodonStart,
                        alignedCodonEnd,
                        prevRefCodon.toLowerCase() + refCodon.toLowerCase() + ">" + prevRefCodon.toLowerCase() + altCodon.toUpperCase() + refCodon.toLowerCase());
            }
        }
        else {

            // Get the capitalization correct in the alt string:
            final String capitalizedAltString = createCapitalizedAlternateAlleleInsertionCodonChangeString(
                    seqComp.getAlignedCodingSequenceAlternateAllele(),
                    seqComp.getAlternateAllele(),
                    seqComp.getCodingSequenceAlleleStart() - seqComp.getAlignedCodingSequenceAlleleStart() -
                            // Need to subtract 1 if we're going on the other strand to line up correctly with
                            // the start of the variant alleles:
                            (seqComp.getStrand() == Strand.POSITIVE ? 0 : 1)
            );

            return String.format(CODON_CHANGE_FORMAT_STRING,
                    alignedCodonStart,
                    alignedCodonEnd,
                    refCodon + ">" + capitalizedAltString);
        }
    }

    private static String getCodonChangeStringForFrameShifts(final SequenceComparison seqComp, final boolean isInsertion, int alignedCodonStart, int alignedCodonEnd, String refCodon) {

        // Requires:
        //     seqComp.getCodingSequenceAlleleStart()
        //     seqComp.getStrand()
        //     seqComp.getTranscriptCodingSequence()

        // Some special deletion adjustments to account for the VCF leading base problem.
        //
        if ( ((seqComp.getCodingSequenceAlleleStart() % AminoAcid.CODON_LENGTH) == 0) ) {
            // If we're an insertion and we start at the begining of a codon, we have to adjust for
            // the preceding base in the input VCF.  This means we have to skip the first or last AminoAcid.CODON_LENGTH bases in the aligned
            // codon, because we have included them to grab the aligned data for the preceding base.
            // Note: This is not the aligned coding sequence position, but the raw coding sequence position:
            if ( isInsertion ) {
                if ( seqComp.getStrand() == Strand.POSITIVE ) {
                    // Skip the first AminoAcid.CODON_LENGTH bases in the aligned codon:
                    alignedCodonStart += AminoAcid.CODON_LENGTH;
                    alignedCodonEnd += AminoAcid.CODON_LENGTH;
                    // Get the next bases in the coding sequence:
                    // TODO: Make sure this won't fail for insertions at the end of a transcript!
                    refCodon = seqComp.getTranscriptCodingSequence()
                            .getBaseString()
                            // Subtract 1 because we're 1-based for genomic coordinates:
                            .substring(alignedCodonStart-1, alignedCodonEnd)
                            .toLowerCase();
                }
            }
            // If we're a deletion and we start at the begining of a codon, we have to adjust for
            // the preceding base in the input VCF.  This means we have to skip the first or last AminoAcid.CODON_LENGTH bases in the aligned
            // codon, because we have included them to grab the aligned data for the preceding base.
            // Note: This is not the aligned coding sequence position, but the raw coding sequence position:
            else {
                if ( seqComp.getStrand() == Strand.POSITIVE ) {
                    // Skip the first AminoAcid.CODON_LENGTH bases in the aligned codon:
                    alignedCodonStart += AminoAcid.CODON_LENGTH;
                    refCodon = refCodon.substring(AminoAcid.CODON_LENGTH);
                }
                else {
                    // Skip the last AminoAcid.CODON_LENGTH bases in the aligned codon:
                    alignedCodonEnd -= AminoAcid.CODON_LENGTH;
                    refCodon = refCodon.substring(0, refCodon.length() - AminoAcid.CODON_LENGTH);
                }
            }
        }

        return String.format(CODON_CHANGE_FORMAT_STRING, alignedCodonStart, alignedCodonEnd, refCodon + "fs");
    }

    private static boolean isStartCodonIndel(final SequenceComparison seqComp, final Locatable startCodon) {
        // Check if the indel is a start codon insertion / deletion.
        // If so, we do not render codon change strings (replicates Oncotator behavior):

        if ( (startCodon != null) &&
             GATKVariantContextUtils.isIndel(seqComp.getReferenceAllele(), seqComp.getAlternateAllele()) ) {

            final SimpleInterval variantInterval = new SimpleInterval(seqComp.getContig(),
                    // Add 1 here because of the leading base for VCFs:
                    seqComp.getAlleleStart() + 1,
                    seqComp.getAlleleStart() + seqComp.getReferenceAllele().length());

            if ( startCodon.overlaps(variantInterval) ) {
                return true;
            }
        }
        return false;
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
     * Gets the next (+ strand) or previous (- strand) complete in-frame codon from the given {@link ReferenceSequence}
     * according to the current codon position and strand.
     * @param referenceSequence The {@link ReferenceSequence} containing the complete coding sequence for the transcript on which the current variant occurs.  Must not be {@code null}.
     * @param currentAlignedCodingSequenceAlleleStart The starting position (1-based, inclusive) of the current codon.  Must be > 0.
     * @param currentAlignedCodingSequenceAlleleStop The ending position (1-based, inclusive) of the current codon.  Must be > 0.
     * @param strand The {@link Strand} on which the current codon resides.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return The next (+ strand) or previous (- strand) codon in frame with the current codon as specified by the given current codon positions.
     */
    private static String getAdjacentReferenceCodon(final ReferenceSequence referenceSequence,
                                                    final int currentAlignedCodingSequenceAlleleStart,
                                                    final int currentAlignedCodingSequenceAlleleStop,
                                                    final Strand strand) {

        Utils.nonNull( referenceSequence );
        ParamUtils.isPositive(currentAlignedCodingSequenceAlleleStart, "Genomic positions must be > 0.");
        ParamUtils.isPositive(currentAlignedCodingSequenceAlleleStop, "Genomic positions must be > 0.");
        assertValidStrand(strand);

        final String nextRefCodon;
        if ( strand == Strand.POSITIVE ) {

            // Add AminoAcid.CODON_LENGTH to get the "next" codon on the - strand:
            final int endex = currentAlignedCodingSequenceAlleleStop + AminoAcid.CODON_LENGTH;

            // Make sure we don't try to get bases after the end of our reference sequence:
            if ( endex >= referenceSequence.getBaseString().length() ) {
                nextRefCodon = "";
            }
            else {
                nextRefCodon = referenceSequence.getBaseString().substring(currentAlignedCodingSequenceAlleleStop, endex);
            }
        }
        else {
            // Make sure we don't try to get bases before the start of our reference sequence:
            if ( currentAlignedCodingSequenceAlleleStart == 1 ) {
                nextRefCodon = "";
            }
            else {
                // Subtract 1 because of 1-inclusive genomic positions
                // Subtract AminoAcid.CODON_LENGTH to get the "next" codon on the - strand:
                nextRefCodon = referenceSequence.getBaseString().substring(currentAlignedCodingSequenceAlleleStart - 1 - AminoAcid.CODON_LENGTH, currentAlignedCodingSequenceAlleleStart - 1);
            }
        }
        return nextRefCodon;
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
     * Creates the string representation of the protein change for the given {@link SequenceComparison}.
     *
     * Specifically, this creates a rendering of the first amino acids to be changed by a variant.
     *
     * Requires that the given {@code seqComp} has the following fields defined with values that are not {@code null}:
     *     proteinChangeInfo
     *     referenceAllele
     *     alternateAllele
     *     alleleStart
     * @param seqComp {@link SequenceComparison} representing the alternate and reference alleles for a DNA sequence.  Must not be {@code null}.
     * @param startCodon {@link Locatable} representing the start codon for the coding region for the variant in {@code seqComp}.  May be {@code null}.
     * @return A {@link String} representing the codon change for the given {@link SequenceComparison}.
     */
    public static String renderProteinChangeString(final SequenceComparison seqComp,
                                                   final Locatable startCodon) {
        Utils.nonNull(seqComp.getProteinChangeInfo());
        Utils.nonNull(seqComp.getReferenceAllele());
        Utils.nonNull(seqComp.getAlternateAllele());

        if ( GATKVariantContextUtils.isIndel(seqComp.getReferenceAllele(), seqComp.getAlternateAllele())) {
            Utils.nonNull(seqComp.getAlleleStart());
        }

        final String refAaSeq        = seqComp.getProteinChangeInfo().getRefAaSeq();
        final String altAaSeq        = seqComp.getProteinChangeInfo().getAltAaSeq();
        final int protChangeStartPos = seqComp.getProteinChangeInfo().getAaStartPos();
        final int protChangeEndPos   = seqComp.getProteinChangeInfo().getAaEndPos();

        // NOTE:
        // It is possible for either amino acid sequence to contain AminoAcid.UNDECODABLE, however that should not
        // actually impact the rendered amino acid sequence other than to put a "?" in the protein sequence in the
        // appropriate place.

        if ( isStartCodonIndel(seqComp, startCodon) ) {
            return "";
        }
        else if ( GATKVariantContextUtils.isFrameshift( seqComp.getReferenceAllele(), seqComp.getAlternateAllele() ) ) {
            return String.format(PROTEIN_CHANGE_FORMAT_STRING, "", refAaSeq, protChangeStartPos, "fs");
        }
        else if (GATKVariantContextUtils.isInsertion( seqComp.getReferenceAllele(), seqComp.getAlternateAllele() )) {
            //  In-Frame Insertion
            if ( protChangeStartPos == protChangeEndPos) {
                return String.format(PROTEIN_CHANGE_FORMAT_STRING, "", refAaSeq, protChangeStartPos, altAaSeq);
            }
            else if ( refAaSeq.isEmpty() ) {
                return String.format(PROTEIN_CHANGE_FORMAT_STRING, "", protChangeStartPos + "_" + protChangeEndPos , "ins", altAaSeq);
            }
            else {
                return String.format(PROTEIN_CHANGE_FORMAT_STRING, protChangeStartPos + "_" + protChangeEndPos, refAaSeq , ">", altAaSeq);
            }
        }
        else if (GATKVariantContextUtils.isDeletion( seqComp.getReferenceAllele(), seqComp.getAlternateAllele() )) {
            //  In-Frame Deletion
            if ( protChangeStartPos != protChangeEndPos ) {
                return String.format(PROTEIN_CHANGE_FORMAT_STRING, protChangeStartPos + "_" + protChangeEndPos, refAaSeq , ">", altAaSeq);
            }
            else {
                return String.format(PROTEIN_CHANGE_FORMAT_STRING, "", refAaSeq , protChangeStartPos, "del");
            }
        }
        else {
            if ( protChangeStartPos != protChangeEndPos ) {
                return String.format(PROTEIN_CHANGE_FORMAT_STRING, protChangeStartPos + "_" + protChangeEndPos, refAaSeq , ">", altAaSeq);
            }
            else {
                return String.format(PROTEIN_CHANGE_FORMAT_STRING, "", refAaSeq , protChangeStartPos, altAaSeq);
            }
        }
    }

    /**
     * Get the coding sequence change string from the given variant information.
     * Also known as cDNA string.
     * This method is assumed to be called only when the variant occurs in a coding region of the genome.
     * @param codingSequenceAlleleStart The start position (1-based, inclusive) of the ref allele in the coding sequence.
     * @param refAllele A {@link String} containing the bases of the reference allele.  MUST BE Reverse Complemented if strand == {@link Strand#NEGATIVE}
     * @param altAllele A {@link String} containing the bases of the alternate allele.  MUST BE Reverse Complemented if strand == {@link Strand#NEGATIVE}
     * @param strand The {@link Strand} on which the gene containing the variant lies.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @param exonStartPosition An {@link Integer} containing the start position (1-based, inclusive) of the variant in its exon.  May be {@code null}.
     * @param exonEndPosition  An {@link Integer} containing the end position (1-based, inclusive) of the variant in its exon.  May be {@code null}.
     * @param alleleStart The start position (1-based, inclusive) of the variant in genomic coordinates (relative to the start of the variant's contig).
     * @return A {@link String} representing the coding sequence change between the given ref and alt alleles.
     */
    public static String getCodingSequenceChangeString( final int codingSequenceAlleleStart,
                                                        final String refAllele,
                                                        final String altAllele,
                                                        final Strand strand,
                                                        final Integer exonStartPosition,
                                                        final Integer exonEndPosition,
                                                        final Integer alleleStart) {
        Utils.nonNull(refAllele);
        Utils.nonNull(altAllele);
        assertValidStrand(strand);

        // Check for ONP:
        if ( GATKVariantContextUtils.isXnp(refAllele, altAllele) ) {
            if (altAllele.length() > 1) {
                return String.format(CDNA_CHANGE_FORMAT_STRING,
                        codingSequenceAlleleStart + "_" + (codingSequenceAlleleStart + refAllele.length() - 1),
                        refAllele,
                        ">",
                        altAllele);
            } else {
                return String.format(CDNA_CHANGE_FORMAT_STRING, codingSequenceAlleleStart, refAllele, ">", altAllele);
            }
        }
        // Check for Insertion:
        else if ( GATKVariantContextUtils.isInsertion(refAllele, altAllele) ) {
            // Must account for strandedness when dealing with insertions:
            // NOTE: Purposefully removing last base in allele to get rid of padding base.
            if ( strand == Strand.NEGATIVE ) {
                return String.format(CDNA_CHANGE_FORMAT_STRING, (codingSequenceAlleleStart-1) + "_" + codingSequenceAlleleStart, "", "ins", altAllele.substring(0,altAllele.length()-1).toUpperCase());
            }
            else {
                return String.format(CDNA_CHANGE_FORMAT_STRING, codingSequenceAlleleStart + "_" + (codingSequenceAlleleStart + 1), "", "ins", altAllele.substring(refAllele.length()).toUpperCase());
            }
        }
        // Must be a Deletion:
        else {
            int start = codingSequenceAlleleStart + altAllele.length();
            int end = codingSequenceAlleleStart + refAllele.length() - 1;

            String deletedBases = refAllele.substring(altAllele.length()).toUpperCase();

            // Must account for strandedness when dealing with deletions:
            if ( strand == Strand.NEGATIVE ) {

                // Because we're backwards stranded we need to account for being on the other side of the reference
                // base, so we shift everything by 1:
                --start;
                --end;

                // We need to adjust our deleted bases here.
                // We know that the last base in the allele is the required preceding reference base.
                deletedBases = refAllele.substring(0, refAllele.length() - 1).toUpperCase();
            }

            // If we have exon information, and we SHOULD, we use it to trim the start/stop coordinates
            // of the cDNA string to the extants of the coding region:
            if ( (exonStartPosition != null) && (exonEndPosition != null) ) {
                final int cdsExonStart = codingSequenceAlleleStart - (alleleStart - exonStartPosition);
                final int cdsExonEnd   = cdsExonStart + (exonEndPosition - exonStartPosition);

                if ( start < cdsExonStart ) {
                    start = cdsExonStart;
                }
                if ( end > cdsExonEnd ) {
                    end = cdsExonEnd;
                }
            }

            // Only add end extent if we have to:
            final String endBoundString = start == end ? "" : "_" + end;

            return String.format(CDNA_CHANGE_FORMAT_STRING, start + endBoundString, "", "del", deletedBases);
        }
    }

    /**
     * Creates an amino acid sequence from a given coding sequence.
     * If the coding sequence is not evenly divisible by {@link AminoAcid#CODON_LENGTH}, the remainder bases will not be included in the coding sequence.
     * If the coding sequence cannot be fully decoded due to IUPAC bases in the ref or alt allele, then this method will
     * include bases representing {@link AminoAcid#UNDECODABLE}.
     * @param codingSequence The coding sequence from which to create an amino acid sequence.  Must not be {@code null}.
     * @param isFrameshift Whether the given {@code codingSequence} was derived from a frameshift mutation.  In this case, no warning will be issued for incorrect sequence length.
     * @param extraLoggingInfo A {@link String} containing extra info for logging purposes.
     * @return A {@link String} containing a sequence of single-letter amino acids.
     */
    static String createAminoAcidSequence(final String codingSequence, final boolean isFrameshift, final String extraLoggingInfo) {
        return createAminoAcidSequenceHelper(codingSequence, isFrameshift, false, extraLoggingInfo);
    }

    /**
     * Creates a Mitochondrial amino acid sequence from a given coding sequence.
     * If the coding sequence is not evenly divisible by 3, the remainder bases will not be included in the coding sequence.
     * If the coding sequence cannot be fully decoded due to IUPAC bases in the ref or alt allele, then this method will
     * include bases representing {@link AminoAcid#UNDECODABLE}.
     * @param codingSequence The coding sequence from which to create an amino acid sequence.  Must not be {@code null}.
     * @param isFrameshift Whether the given {@code codingSequence} was derived from a frameshift mutation.  In this case, no warning will be issued for incorrect sequence length.
     * @param extraLoggingInfo A {@link String} containing extra info for logging purposes.
     * @return A {@link String} containing a sequence of single-letter amino acids.
     */
    static String createMitochondrialAminoAcidSequence(final String codingSequence, final boolean isFrameshift, final String extraLoggingInfo) {
        return createAminoAcidSequenceHelper(codingSequence, isFrameshift, true, extraLoggingInfo);
    }

    private static String createAminoAcidSequenceHelper(final String codingSequence, final boolean isFrameshift, final boolean isMitochondria, final String extraLoggingInfo) {
        Utils.nonNull(codingSequence);

        final StringBuilder sb = new StringBuilder();

        // Ensure that we don't have remainder bases:
        int maxIndex = codingSequence.length();
        if ( maxIndex % AminoAcid.CODON_LENGTH != 0 ) {
            maxIndex = (int)Math.floor(maxIndex / AminoAcid.CODON_LENGTH) * AminoAcid.CODON_LENGTH;
            if ( !isFrameshift ) {
                logger.warn("createAminoAcidSequence given a coding sequence of length not divisible by " + AminoAcid.CODON_LENGTH + ".  Dropping bases from the end: " + (codingSequence.length() % AminoAcid.CODON_LENGTH) + (extraLoggingInfo.isEmpty() ? "" : " " + extraLoggingInfo));
            }
        }

        // Call the correct method based on whether or not we're looking up mitochondria:
        // NOTE: Based on conversations with a Mitochondria expert (Sarah Calvo), the alternate start codons should not
        // be taken into account here, so we don't have to use the overload with the booleans.
        final Function<String, AminoAcid> aminoAcidLookupFunction = ( isMitochondria ? FuncotatorUtils::getMitochondrialAminoAcidByCodon : FuncotatorUtils::getEukaryoticAminoAcidByCodon );

        for ( int i = 0; i < maxIndex; i += AminoAcid.CODON_LENGTH ) {
            final AminoAcid aa = aminoAcidLookupFunction.apply(codingSequence.substring(i, i+3));
            if ( aa == null ) {
                logger.warn("Unexpected codon sequence: " + codingSequence.substring(i, i+AminoAcid.CODON_LENGTH) + "   -   Cannot create protein prediction.  This may be due to \"N\" (or other IUPAC) bases in the Reference or the current variant's alleles.");

                // Because we couldn't decode the amino acid we must append the UNDECODABLE AminoAcid sequence:
                sb.append(AminoAcid.UNDECODABLE.getLetter());
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

        // NOTE: This check appears to be redundant, but in actuality, it is required.
        //       Because we reconstruct the coding sequence allele separately from the reference allele, we need to check this
        //       again to make sure we have the right alleles given our input.
        if ( !expectedReferenceSequence.equals(refAllele.getBaseString()) ) {
            // Oh noes!
            // Ref allele is different from reference sequence!
            // Oh well, we should use the reference we were given anyways...
            final String substitutedAlignedSeq = getAlternateSequence(new StrandCorrectedReferenceBases(codingSequence, strand), refAlleleStart, refAllele, refAllele, strand);

            // We use the positive strand here because we have already reverse complemented the sequence in the call
            // above.
            final String substitutedAlignedAlleleSeq = getAlignedAlleleSequence(substitutedAlignedSeq, alignedAlleleStart, alignedAlleleStop, Strand.POSITIVE);

            // Warn the user!
            logger.warn("Reference allele differs from reference coding sequence (strand: " + strand + ") (Allele " + expectedReferenceSequence + " != " + refAllele.getBaseString() + " coding sequence)!  Substituting given allele for sequence code (" + alignedAlleleSeq + "->" + substitutedAlignedAlleleSeq + ")");

            // Set up our return value:
            alignedAlleleSeq = substitutedAlignedAlleleSeq;
        }

        return alignedAlleleSeq;
    }

    /**
     * Get the aligned coding sequence for the given reference allele.
     * NOTE: alignedRefAlleleStart must be less than or equal to codingSequenceRefAlleleStart.
     * @param referenceSnippet {@link StrandCorrectedReferenceBases} containing a short excerpt of the reference sequence around the given reference allele.  Must not be {@code null}.
     * @param referencePadding Number of bases in {@code referenceSnippet} before the reference allele starts.  This padding exists at the end as well (plus some other bases to account for the length of the alternate allele if it is longer than the reference).  Must be >= 0.
     * @param refAllele The reference {@link Allele}.  Must not be {@code null}.
     * @param altAllele The alternate {@link Allele}.  Must not be {@code null}.
     * @param codingSequenceRefAlleleStart The position (1-based, inclusive) in the coding sequence where the {@code refAllele} starts.
     * @param alignedRefAlleleStart The in-frame position (1-based, inclusive) of the first base of the codon containing the reference allele.
     * @param strand The {@link Strand} on which the variant occurs.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @param variantGenomicPositionForLogging The {@link Locatable} representing the position of the variant in genomic coordinates.  Used for logging purposes only.
     * @return A {@link String} of in-frame codons that contain the entire reference allele.
     */
    public static String getAlignedRefAllele(final StrandCorrectedReferenceBases referenceSnippet,
                                             final int referencePadding,
                                             final Allele refAllele,
                                             final Allele altAllele,
                                             final int codingSequenceRefAlleleStart,
                                             final int alignedRefAlleleStart,
                                             final Strand strand,
                                             final Locatable variantGenomicPositionForLogging) {

        Utils.nonNull(referenceSnippet);
        Utils.nonNull(refAllele);
        Utils.nonNull(altAllele);
        ParamUtils.isPositiveOrZero( referencePadding, "Padding must be >= 0." );
        Utils.nonNull(variantGenomicPositionForLogging);
        assertValidStrand(strand);
        Utils.validate( alignedRefAlleleStart <= codingSequenceRefAlleleStart, "The alignedRefAlleleStart must be less than or equal to codingSequenceRefAlleleStart!" );

        // Get the length of the variant.
        // Since our referenceSnippet is <referencePadding> bases on either side of the variant
        // we need to know how long the variant is (so we can adjust for - strand variants)
        final int referenceLength = refAllele.length();

        // Compute the number of bases that we must prepend to this one to get to the start of a codon:
        final int extraBasesNeededForCodonAlignment = (codingSequenceRefAlleleStart - alignedRefAlleleStart);

        // We must add bases to the start to adjust for indels because of the required preceding base in VCF format:
        final int indelAdjustment = GATKVariantContextUtils.isIndel(refAllele, altAllele) ? 1 : 0;

        // We need to adjust coordinates for deletions on the - strand in the reference snippet:
        final int negativeDeletionAdjustment = (altAllele.length() < refAllele.length()) && (strand == Strand.NEGATIVE) ? 1 : 0;

        // Get the index in our reference snippet of the aligned codon start:
        int refStartAlignedIndex;
        if ( strand == Strand.NEGATIVE ) {
            // Variants on the - strand have to be handled as a special case because of how
            // the referenceSnippet is constructed:
            refStartAlignedIndex = referencePadding - extraBasesNeededForCodonAlignment;
        }
        else {
            refStartAlignedIndex = referencePadding - extraBasesNeededForCodonAlignment - indelAdjustment;
        }

        // TODO: This should probably be an error condition:
        if ( refStartAlignedIndex < 0 ) {
            refStartAlignedIndex = 0;
        }

        // Round to the nearest multiple of AminoAcid.CODON_LENGTH to get the end position.
        // Note this is the position of the first base NOT in the REF allele
        // (because it's used for substring coordinates).
        int refEndPosExclusive = refStartAlignedIndex + (int)(Math.ceil((extraBasesNeededForCodonAlignment + refAllele.length()) / ((double)AminoAcid.CODON_LENGTH)) * AminoAcid.CODON_LENGTH);

        // Create the aligned reference:
        String alignedReferenceAllele = referenceSnippet.getBaseString().substring(refStartAlignedIndex, refEndPosExclusive);

        final String computedReferenceAlleleWithSubstitution = alignedReferenceAllele.substring(extraBasesNeededForCodonAlignment, extraBasesNeededForCodonAlignment + refAllele.length());

        // Make sure our reference is what we expect it to be with the ref allele:
        if ( !computedReferenceAlleleWithSubstitution.equals(refAllele.getBaseString()) ) {
            // Oh noes!
            // Ref allele is different from reference sequence!

            // Oh well, we should use the reference we were given anyways...
            final String substitutedReferenceSnippet = getAlternateSequence(referenceSnippet, referencePadding + 1, refAllele, refAllele, strand);
            refEndPosExclusive = refStartAlignedIndex + (int)(Math.ceil((extraBasesNeededForCodonAlignment + refAllele.length()) / ((double)AminoAcid.CODON_LENGTH)) * AminoAcid.CODON_LENGTH);

            final String substitutedAlignedAlleleSeq = substitutedReferenceSnippet.substring(refStartAlignedIndex, refEndPosExclusive);

            // Warn the user!
            final String positionString = '[' + variantGenomicPositionForLogging.getContig() + ":" + variantGenomicPositionForLogging.getStart() + ']';
            logger.warn("Reference allele is different than the reference coding sequence (strand: " + strand + ", alt = " + altAllele.getBaseString() + ", ref " + refAllele.getBaseString() + " != " + computedReferenceAlleleWithSubstitution + " reference coding seq) @" + positionString + "!  Substituting given allele for sequence code (" + alignedReferenceAllele + "->" + substitutedAlignedAlleleSeq + ")");

            // Set up our return value:
            alignedReferenceAllele = substitutedAlignedAlleleSeq;
        }

        return alignedReferenceAllele;
    }

    /**
     * Get the bases around the given ref allele in the correct direction of the strand for this variant.
     * The number of bases before and after the variant is specified by {@code referenceWindow}.
     * The result will be trimmed down by 1 base in the event that the variant is an indel in order to account for
     * the required preceding base in VCF format.
     *
     * ASSUMES: that the given {@link ReferenceContext} is already centered on the variant location.
     *
     * @param refAllele The reference {@link Allele} for the variant.  Used
     * @param altAllele The alternate {@link Allele} for the variant.
     * @param reference The {@link ReferenceContext} for the variant, with the current window centered on the variant with no padding around it.
     * @param strand The {@link Strand} on which the variant occurs.
     * @param referenceWindow The number of bases on either side of the variant to add to the resulting string.
     * @return A {@link StrandCorrectedReferenceBases} of bases of length either {@code referenceWindow} * 2 + |ref allele| OR {@code referenceWindow} * 2 + |ref allele| - 1 , corrected for strandedness.
     */
    public static StrandCorrectedReferenceBases createReferenceSnippet(final Allele refAllele, final Allele altAllele, final ReferenceContext reference, final Strand strand, final int referenceWindow ) {

        // Make sure our window in the reference includes only our variant:
        Utils.validate(
                (reference.numWindowLeadingBases() == reference.numWindowTrailingBases()) && (reference.numWindowLeadingBases() == 0),
                "Reference must have no extra bases around the variant.  Found: " + reference.numWindowLeadingBases() + " before / + "  + reference.numWindowTrailingBases() + " + after"
        );

        // We must add bases to the start to adjust for indels because of the required preceding base in VCF format:
        final int indelStartBaseAdjustment = GATKVariantContextUtils.isIndel(refAllele, altAllele) ? 1 : 0;

        final int start = reference.getWindow().getStart() - referenceWindow + indelStartBaseAdjustment;
        final int end   = reference.getWindow().getEnd() + referenceWindow;

        // Calculate the interval from which to get the reference:
        final SimpleInterval refBasesInterval = new SimpleInterval(reference.getWindow().getContig(), start, end);

        // Get the reference bases for this interval.
        final byte[] referenceBases = reference.getBases(refBasesInterval);

        // Get the bases in the correct direction:
        if ( strand == Strand.POSITIVE ) {
            return new StrandCorrectedReferenceBases(referenceBases, strand);
        }
        else {
            return new StrandCorrectedReferenceBases(ReadUtils.getBasesReverseComplement(referenceBases), strand);
        }
    }

    /**
     * Get a string of bases around a variant (specified by reference and alternate alleles), including the reference allele itself.
     * ASSUMES: that the given {@link ReferenceContext} is already centered on the variant location.
     * @param refAllele The reference {@link Allele} for the variant in question.  If on {@link Strand#NEGATIVE}, must have already been reverse complemented.  Must not be {@code null}.
     * @param referenceContext The {@link ReferenceContext} centered around the variant in question.  Must not be {@code null}.
     * @param strand The {@link Strand} on which the variant in question lives.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @param referenceWindowInBases The number of bases to the left and right of the variant to return.  Must be > 0.
     * @return A {@link StrandCorrectedReferenceBases} containing {@code referenceWindowInBases} bases to either side of the specified refAllele.
     */
    public static StrandCorrectedReferenceBases getBasesInWindowAroundReferenceAllele( final Allele refAllele,
                                                                final ReferenceContext referenceContext,
                                                                final Strand strand,
                                                                final int referenceWindowInBases) {
        // Call into the other method that we have to do this sort of thing.
        // We use the ref allele as both ref and alt so that we can avoid adding indel bases:
        return createReferenceSnippet( refAllele, refAllele, referenceContext, strand, referenceWindowInBases );
    }

    /**
     * Get the Protein change start position (1-based, inclusive) given the aligned position of the coding sequence.
     * @param alignedCodingSequenceAlleleStart Position (1-based, inclusive) of the start of the allele in the coding sequence.  Must not be {@code null}.  Must be > 0.
     * @return The position (1-based, inclusive) of the protein change in the amino acid sequence.
     */
    public static int getProteinChangePosition(final Integer alignedCodingSequenceAlleleStart) {

        Utils.nonNull(alignedCodingSequenceAlleleStart);
        ParamUtils.isPositive( alignedCodingSequenceAlleleStart, "Genome positions must be > 0." );

        return ((alignedCodingSequenceAlleleStart-1) / AminoAcid.CODON_LENGTH) + 1; // Add 1 because we're 1-based.
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
     * @param altAllele Alternate Allele.  Used for both content and length of the alternate allele.  If on the negative strand, assumes that the bases are already reverse-complemented, and that the leading reference base is the last base in the allele.  Must not be {@code null}.
     * @param strand The {@link Strand} on which the variant occurs.  Must not be {@code null}.  Must not be {@link Strand#NONE}.
     * @return The coding sequence that includes the given alternate allele in place of the given reference allele.
     */
    public static String getAlternateSequence(final StrandCorrectedReferenceBases referenceSequence,
                                              final int alleleStartPos,
                                              final Allele refAllele,
                                              final Allele altAllele,
                                              final Strand strand) {

        Utils.nonNull(referenceSequence);
        Utils.nonNull(refAllele);
        Utils.nonNull(altAllele);

        ParamUtils.isPositive( alleleStartPos, "Genome positions must be > 0." );

        // We have to subtract 1 here because we need to account for the 1-based indexing of
        // the start and end of the coding region:
        final int alleleIndex = Math.abs(alleleStartPos - 1);

        return referenceSequence.getBaseString().substring(0, alleleIndex) +
                altAllele.getBaseString() +
                referenceSequence.getBaseString().substring(alleleIndex + refAllele.length());
    }

    /**
     * Get the position (1-based, inclusive) of the given {@link VariantContext} start relative to the transcript it appears in.
     * The transcript is specified by {@code sortedTranscriptExonList}.
     * @param variant The {@link VariantContext} of which to find the start position in the given transcript (must not be {@code null}).
     * @param exons {@link List} of {@link Locatable}s representing the exons in the transcript in which the given {@code variant} occurs.
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
        final int numExonsBeforeLast = sortedFilteredExons.size() - 1;
        for ( int i = 0; i < numExonsBeforeLast; ++i ) {
            final Locatable exon = sortedFilteredExons.get(i);

            // Add 1 to position to account for inclusive values:
            position += exon.getEnd() - exon.getStart() + 1;
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

        if ( B37_SEQUENCE_DICTIONARY == null ) {
            B37_SEQUENCE_DICTIONARY = initializeB37SequenceDict();
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
        return B37_To_HG19_CONTIG_NAME_MAP.getOrDefault(b37Contig, b37Contig);
    }

    /**
     * Converts a given HG19 style contig name to the equivalent in B37.
     * @param hg19Contig The contig name from the HG19 Human Genome reference to convert to the equivalent contig from the B37 Human Genome reference.
     * @return The B37 equivalent of the given HG19 contig name, if such an equivalent name exists.  If no equivalent name exists, returns the given {@code hg19Contig}.
     */
    public static String convertHG19ContigToB37Contig( final String hg19Contig ) {
        return HG19_TO_B37_CONTIG_NAME_MAP.getOrDefault(hg19Contig, hg19Contig);
    }

    /**
     * Create a cDNA string for an intronic variant.
     * The cDNA string contains information about the coding position of a variant and the alleles involved.
     * If the transcript in which the variant occurs contains at least 1 exon, the string will be non-empty and of
     * the form:
     *     c.e[EXON NUMBER][+|-][BASES FROM EXON][REF ALLELE]>[ALT ALLELE]
     * Concretely:
     *     c.e2-1A>G
     * Where:
     *      2 = the number of the exon to which the given variant start is closest
     *     -1 = number of bases away from the exon (1 before)
     *      A = Reference allele
     *      G = Alternate allele
     * @param variantStart The start position (1-based, inclusive) in genomic coordinates of the variant.
     * @param exonList The {@link List<Locatable>} representing the exons in the transcript in which this variant occurs.  Must not be {@code null}.
     * @param strandCorrectedRefAllele A {@link String} containing the bases of the reference allele, which are correct for strandedness (i.e. if on a transcript on the {@link Strand#NEGATIVE} strand, the string has already been reverse-complemented).  Must not be {@code null}.
     * @param strandCorrectedAltAllele A {@link String} containing the bases of the alternate allele, which are correct for strandedness (i.e. if on a transcript on the {@link Strand#NEGATIVE} strand, the string has already been reverse-complemented).  Must not be {@code null}.
     * @return A {@link String} representing the cDNA change for the given data.  Will be empty if the given {@code exonList} is empty.
     */
    public static String createIntronicCDnaString(final int variantStart,
                                                  final List<? extends Locatable> exonList,
                                                  final String strandCorrectedRefAllele,
                                                  final String strandCorrectedAltAllele) {

        Utils.nonNull(exonList);
        Utils.nonNull(strandCorrectedRefAllele);
        Utils.nonNull(strandCorrectedAltAllele);

        // Get the exon that is closest to our variant:
        final int exonIndex = getClosestExonIndex(variantStart, exonList);

        if ( exonIndex != -1 ) {
            final Locatable closestExon = exonList.get(exonIndex);

            final int startDiff = variantStart - closestExon.getStart();
            final int endDiff =  variantStart - closestExon.getEnd();

            // Get the offset from our start:
            final int exonOffset;
            if ( Math.abs(startDiff) <= Math.abs(endDiff) ) {
                exonOffset = startDiff;
            }
            else {
                exonOffset = endDiff;
            }

            // Get the cDNA string itself:
            return "c.e" + (exonIndex+1) + (exonOffset < 0 ? "-" : "+") + Math.abs(exonOffset) + strandCorrectedRefAllele + ">" + strandCorrectedAltAllele;
        }
        else {
            return "NA";
        }
    }


    /**
     * Get the index of the exon that is closest to the given start position of a variant.
     * Checks both before and after the start position to get the closest exon such that it may occur before the
     * variant, within the variant, or after the variant.
     * If there are no exons in the transcript, will return {@code null}.
     * @param variantStartPos Start position (1-based, inclusive) in genomic coordinates of a variant.
     * @param exonList The {@link List<Locatable>} representing the exons in the transcript in which this variant occurs.  Must not be {@code null}.
     * @return The index into the given {@code exonList} corresponding to the entry which is the fewest bases away from the given variant positions.  If {@code exonList} is empty, will return {@code -1}.
     */
     public static int getClosestExonIndex( final int variantStartPos,
                                             final List<? extends Locatable> exonList) {
        Utils.nonNull(exonList);

        int exonIndex = -1;
        int distFromVariant = Integer.MAX_VALUE;

        for ( int i = 0; i < exonList.size() ; ++i ) {
            for ( final int exonPos : Arrays.asList(exonList.get(i).getStart(), exonList.get(i).getEnd()) ) {
                final int dist = Math.abs(variantStartPos - exonPos);
                if ( dist < distFromVariant ) {
                    exonIndex = i;
                    distFromVariant = dist;
                }
            }
        }

        return exonIndex;
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
     * @return A {@link SimpleInterval} containing the extents of the exon that overlaps the given altAllele.  {@code null} if no overlap occurs.
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
    private static synchronized final SAMSequenceDictionary initializeB37SequenceDict() {

        if ( B37_SEQUENCE_DICTIONARY == null ) {
            try {
                final File b37SeqDictFile = Resource.getResourceContentsAsFile("org/broadinstitute/hellbender/tools/funcotator/Homo_sapiens_assembly19.dict");

                return ReferenceUtils.loadFastaDictionary(b37SeqDictFile);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Unable to load b37 dict from jar resources due to IO Exception!", ex);
            }
        }
        else {
            return B37_SEQUENCE_DICTIONARY;
        }
    }

    /**
     * Determines whether the given {@code funcotation} has a transcript ID that is in the given {@code acceptableTranscripts}.
     * Ignores transcript version numbers.
     * @param funcotation The {@link GencodeFuncotation} to check against the set of {@code acceptableTranscripts}.
     * @param acceptableTranscripts The {@link Set} of transcript IDs that are OK to keep.
     * @return {@code true} if funcotation.annotationTranscript is in {@code acceptableTranscripts} (ignoring transcript version); {@code false} otherwise.
     */
    public static boolean isFuncotationInTranscriptList( final GencodeFuncotation funcotation,
                                                  final Set<String> acceptableTranscripts ) {
        if ( funcotation.getAnnotationTranscript() != null ) {
            final List<String> acceptableTranscriptsWithoutVersionNumbers = acceptableTranscripts.stream().map(FuncotatorUtils::getTranscriptIdWithoutVersionNumber).collect(Collectors.toList());
            return acceptableTranscriptsWithoutVersionNumbers.contains( getTranscriptIdWithoutVersionNumber(funcotation.getAnnotationTranscript()) );
        }
        else {
            return false;
        }
    }

    /**
     * Removes the transcript ID version number from the given transcript ID (if it exists).
     * @param transcriptId The transcript from which to remove the version number.
     * @return The {@link String} corresponding to the given {@code transcriptId} without a version number.
     */
    public static String getTranscriptIdWithoutVersionNumber( final String transcriptId ) {
        return transcriptId.replaceAll("\\.\\d+$", "");
    }

    /**
     * Create a configuration instance from a properties file, where key names are ordered.
     *
     * @param configFile Must be a readable file.
     * @return Configuration instance. Never {@code null}
     */
    public static Configuration retrieveConfiguration(final File configFile) {
        IOUtils.assertFileIsReadable(configFile.toPath());

        // Use Apache Commons configuration since it will preserve the order of the keys in the config file.
        //  Properties will not preserve the ordering, since it is backed by a HashSet.
        try {
            return new Configurations().properties(configFile);
        } catch (final ConfigurationException ce) {
            throw new UserException.BadInput("Unable to read from XSV config file: " + configFile.toString(), ce);
        }
    }

    /** Determine whether a given variant context represents a segment.
     *
     * Dev note: this is done by examining the length and the alt allele of the segment.
     *
     * @param vc Never {@code null}
     * @param minSizeForSegment Minimum size for a segment to be valid.
     * @return Boolean whether the given variant context could represent a copy number segment.
     */
    public static boolean isSegmentVariantContext(final VariantContext vc, final int minSizeForSegment) {
        Utils.nonNull(vc);
        final List<String> ACCEPTABLE_ALT_ALLELES = Stream.concat(
                Stream.of(SimpleSVType.SupportedType.values())
                        .map(s -> SimpleSVType.createBracketedSymbAlleleString(s.toString())),
                Stream.of(AnnotatedIntervalToSegmentVariantContextConverter.COPY_NEUTRAL_ALLELE.getDisplayString(), Allele.UNSPECIFIED_ALTERNATE_ALLELE_STRING)
        ).collect(Collectors.toList());

        boolean acceptableAlternateAllele = false;
        for (final Allele a: vc.getAlternateAlleles()) {
            final String representation = a.getDisplayString();
            if (ACCEPTABLE_ALT_ALLELES.contains(representation)) {
                acceptableAlternateAllele = true;
            }
        }

        return acceptableAlternateAllele && (VariantContextUtils.getSize(vc) > minSizeForSegment);
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

    /**
     * @param funcotationHeaderDescription The raw description of the funcotation info field.  Never {@code null}
     * @return Array of the keys, in proper order.  Never {@code null}
     */
    public static String[] extractFuncotatorKeysFromHeaderDescription(final String funcotationHeaderDescription) {
        Utils.nonNull(funcotationHeaderDescription);

        final String[] descriptionSplit = StringUtils.splitByWholeSeparatorPreserveAllTokens(funcotationHeaderDescription,
                VcfOutputRenderer.DESCRIPTION_PREAMBLE_DELIMITER);
        return StringUtils.splitByWholeSeparatorPreserveAllTokens(descriptionSplit[1], VcfOutputRenderer.FIELD_DELIMITER);
    }

    /**
     * Make sure that an individual funcotation field (i.e. single value of a funcotation) is sanitized for VCF consumption.
     * Particularly, make sure that it does not allow special characters that would interfere with VCF parsing.
     * @param individualFuncotationField  value from a funcotation. Never {@code null}
     * @return input string with special characters replaced by _%HEX%_ where HEX is the 2 digit ascii hex code.  Never {@code null}
     */
    static String sanitizeFuncotationFieldForVcf(final String individualFuncotationField) {
        Utils.nonNull(individualFuncotationField);

        // List of letters to encode:
        final List<String> badLetters = Arrays.asList(",", ";", "=", "\t", VcfOutputRenderer.HEADER_LISTED_FIELD_DELIMITER, " ", "\n", VcfOutputRenderer.ALL_TRANSCRIPT_DELIMITER);

        // Encoded version:
        final List<String> cleanLetters = new ArrayList<>();
        for(String s: badLetters) {
            final String hex;
            if (s.getBytes(StandardCharsets.US_ASCII)[0] < 16) {
                hex = "_%0"+Integer.toHexString(s.getBytes(StandardCharsets.US_ASCII)[0] & 0xFF).toUpperCase()+ "_";
            }
            else {
                hex = "_%" +Integer.toHexString(s.getBytes(StandardCharsets.US_ASCII)[0] & 0xFF).toUpperCase() + "_";
            }
            cleanLetters.add(hex);
        }

        // Now replace them:
        return StringUtils.replaceEach(individualFuncotationField, badLetters.toArray(new String[]{}), cleanLetters.toArray(new String[]{}));
    }

    /**
     * Make sure that an individual funcotation field (i.e. single value of a funcotation) is sanitized for MAF consumption.
     * Particularly, make sure that it does not allow special characters that would interfere with MAF parsing.
     * @param individualFuncotationField  value from a funcotation. Never {@code null}
     * @return input string with special characters replaced by _%HEX%_ where HEX is the 2 digit ascii hex code.  Never {@code null}
     */
    public static String sanitizeFuncotationFieldForMaf(final String individualFuncotationField) {
        Utils.nonNull(individualFuncotationField);
        return StringUtils.replaceEach(individualFuncotationField, new String[]{"\t", "\n"}, new String[]{"_%09_", "_%0A_"});
    }

    /**
     * Create a mapping for a single variant.  The mapping are the variant allele(s) to a FuncotationMap {@link FuncotationMap}
     * This is lossy, since the attribute cannot possibly store the type of funcotation.
     *
     * @param funcotationHeaderKeys {@link FuncotatorUtils#extractFuncotatorKeysFromHeaderDescription(String)} for
     *                                  getting this parameter from a VCFHeader entry.  Never {@code null}
     * @param v the variant to use in creating the map.  Never {@code null}
     * @param transcriptIdFuncotationName The field name to use for determining the transcript ID.  Use {@link FuncotationMap#NO_TRANSCRIPT_AVAILABLE_KEY} if unknown.
     *                            If not in the funcotation keys, then the funcotation maps will be created with one transcript ID,
     *                            {@link FuncotationMap#NO_TRANSCRIPT_AVAILABLE_KEY}.  Never {@code null}
     * @param dummyDatasourceName Datasource name to use for the funcotations coming from the VCF.  Note that the original datasource names are impossible to reconstruct.
     *                            Never {@code null}
     * @return Never {@code null}
     */
    public static Map<Allele, FuncotationMap> createAlleleToFuncotationMapFromFuncotationVcfAttribute(final String[] funcotationHeaderKeys,
                                                                                                      final VariantContext v,
                                                                                                      final String transcriptIdFuncotationName,
                                                                                                      final String dummyDatasourceName) {
        Utils.nonNull(funcotationHeaderKeys);
        Utils.nonNull(v);
        Utils.nonNull(transcriptIdFuncotationName);
        Utils.nonNull(dummyDatasourceName);
        final List<String> funcotationPerAllele = v.getAttributeAsList(VcfOutputRenderer.FUNCOTATOR_VCF_FIELD_NAME).stream().map(Object::toString).collect(Collectors.toList());
        if (v.getAlternateAlleles().size() != funcotationPerAllele.size()) {
            throw new GATKException.ShouldNeverReachHereException("Could not parse FUNCOTATION field properly.");
        }

        return IntStream.range(0, v.getAlternateAlleles().size()).boxed()
                .collect(Collectors
                        .toMap(v::getAlternateAllele, i -> FuncotationMap.createAsAllTableFuncotationsFromVcf(transcriptIdFuncotationName,
                                funcotationHeaderKeys, funcotationPerAllele.get(i), v.getAlternateAllele(i), dummyDatasourceName)));
    }

    /**
     * @param f Never {@code null}
     * @return whether this is an instance of {@link GencodeFuncotation}
     */
    public static boolean isGencodeFuncotation(final Funcotation f) {
        Utils.nonNull(f);
        return f instanceof GencodeFuncotation;
    }

    /**
     * @param funcotations  Never {@code null}
     * @return whether any funcotations in the input are an instance of {@link GencodeFuncotation}
     */
    public static boolean areAnyGencodeFuncotation(final List<Funcotation> funcotations) {
        return funcotations.stream().anyMatch(FuncotatorUtils::isGencodeFuncotation);
    }

    /**
     * Create funcotations (one for each alt allele) corresponding to the given variant context.
     *
     * Assumes that the fields in the variant context are named exactly the same as what is in the metadata, though the
     *  metadata may have additional fields.  The metadata must include all variant attributes.
     *
     * @param vc The variant context to derive funcotations.  Never {@code null}
     * @param metadata Existing metadata that must be a superset of the variant context info field attributes.  Never {@code null}
     * @param datasourceName Name to use as the datasource in the funcotations.  Never {@code null}
     * @return A list of funcotations based on the variant context (INFO) attributes.  Never empty, unless the metadata has no fields.  Never {@code null}
     */
    public static List<Funcotation> createFuncotations(final VariantContext vc, final FuncotationMetadata metadata, final String datasourceName) {

        Utils.nonNull(vc);
        Utils.nonNull(metadata);
        Utils.nonNull(datasourceName);

        final List<String> allFields = metadata.retrieveAllHeaderInfo().stream().map(VCFInfoHeaderLine::getID).collect(Collectors.toList());

        final Set<String> attributesNotInMetadata = vc.getAttributes().keySet().stream().filter(k -> !allFields.contains(k)).collect(Collectors.toSet());
        if (attributesNotInMetadata.size() != 0) {
            throw new UserException.MalformedFile("Not all attributes in the variant context appear in the metadata: " + attributesNotInMetadata.stream().collect(Collectors.joining(", ")) + " .... Please add these attributes to the input metadata (e.g. VCF Header).");
        }

        return createFuncotationsFromMetadata(vc, metadata, datasourceName);
    }

    /**
     * Use the given metadata to create funcotations from variant context attributes (and alt alleles)
     * @param vc Never {@code null}
     * @param metadata Fields that should be present in the funcotations.  Can be a superset of the fields in the
     *                 funcotations.  Never {@code null}
     * @param datasourceName Name to appear in all funcotations.  Never {@code null}
     * @return Instances of {@link Funcotation} for each field in the metadata x alternate allele in the variant context.
     * If a field is not present in the variant context attributes, the field will ave value empty string ("") in all output
     * funcotations.  Fields will be the same names and values for each alternate allele in the funcotations.
     */
    static List<Funcotation> createFuncotationsFromMetadata(final VariantContext vc, final FuncotationMetadata metadata, final String datasourceName) {

        Utils.nonNull(vc);
        Utils.nonNull(metadata);
        Utils.nonNull(datasourceName);

        final List<String> fields = metadata.retrieveAllHeaderInfo().stream().map(VCFInfoHeaderLine::getID).collect(Collectors.toList());
        final List<Funcotation> result = new ArrayList<>();
        for (final Allele allele: vc.getAlternateAlleles()) {

            // We must have fields for everything in the metadata.
            final List<String> funcotationFieldValues = new ArrayList<>();
            for (final String funcotationFieldName : fields) {
                funcotationFieldValues.add(vc.getAttributeAsString(funcotationFieldName, ""));
            }

            result.add(TableFuncotation.create(fields, funcotationFieldValues, allele, datasourceName, metadata));
        }

        return result;
    }

    /**
     * @param funcotation Funcotation to render for a VCF.  Never {@code null}
     * @param includedFields List of fields to include.  Any that match fields in the funcotation will be rendered.
     *                       Never {@code null}
     * @return string with the VCF representation of a {@link VcfOutputRenderer#FUNCOTATOR_VCF_FIELD_NAME} for the given
     * Funcotation.  Never {@code null}, but empty string is possible.
     */
    public static String renderSanitizedFuncotationForVcf(final Funcotation funcotation, final List<String> includedFields) {
        Utils.nonNull(funcotation);
        Utils.nonNull(includedFields);
        if (includedFields.size() == 0) {
            return "";
        }
        return funcotation.getFieldNames().stream()
                .filter(includedFields::contains)
                .map(field -> FuncotatorUtils.sanitizeFuncotationFieldForVcf(funcotation.getField(field)))
                .collect(Collectors.joining(VcfOutputRenderer.FIELD_DELIMITER));
    }

    /**
     * A type to keep track of different specific genuses.
     */
    public enum Genus {
        /** Cows */
        BOS,
        /** Humans */
        HOMO,
        /** Mice */
        MUS,
        /** Pheasant */
        CORTURNIX,
        /** Chicken */
        GALLUS,
        /** Unspecified genus / genus not in this list. */
        UNSPECIFIED;
    }

    /**
     * Create a linked hash map from two lists with corresponding values.
     * @param keys Never {@code null}.  Length must be the same as values.
     * @param values Never {@code null}.  Length must be the same as keys.
     * @param <T> Type for keys
     * @param <U> Type for values
     * @return Never {@code null}
     */
    public static <T,U> LinkedHashMap<T,U> createLinkedHashMapFromLists(final List<T> keys, final List<U> values) {
        Utils.nonNull(keys);
        Utils.nonNull(values);
        Utils.validateArg(keys.size() == values.size(), "Keys and values were not the same length.");
        return IntStream.range(0, keys.size()).boxed().collect(Collectors.toMap(keys::get,
                values::get, (x1, x2) -> {
                    throw new IllegalArgumentException("Should not be able to have duplicate field names: " + x1 + " == " + x2);
                }, LinkedHashMap::new));
    }

}