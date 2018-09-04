package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.reflect.Array;
import java.util.List;
import java.util.Objects;
import java.util.OptionalInt;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Miscellaneous utilities to work with Variant data structures.
 */
public class GATKProtectedVariantContextUtils {

    /**
     * Composes the double array from a genotype annotation.
     *
     * @param variantContext the target variant-context.
     * @param key the name of the attribute containing the double array.
     * @param defaultValue the double array to return in case there is no such an annotation.
     * @param missingValue value to use to fill up positions with a missing value (e.g. '.').
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static double[] getAttributeAsDoubleArray(final VariantContext variantContext, final String key,
                                                     final Supplier<double[]> defaultValue, final double missingValue) {
        Utils.nonNull(variantContext);
        return attributeValueToDoubleArray(variantContext.getAttribute(key), key, defaultValue, missingValue);
    }

    /**
     * Composes the double array from a genotype annotation.
     *
     * @param genotype the target variant-context.
     * @param key the name of the attribute containing the double array.
     * @param defaultValue the double array to return in case there is no such an annotation.
     * @param missingValue value to use to fill up positions with a missing value (e.g. '.').
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static double[] getAttributeAsDoubleArray(final Genotype genotype, final String key,
                                                     final Supplier<double[]> defaultValue, final double missingValue) {
        Utils.nonNull(genotype);
        return attributeValueToDoubleArray(genotype.getExtendedAttribute(key), key, defaultValue, missingValue);
    }

    /**
     * Composes the double array from a genotype annotation.
     *
     * @param genotype the target variant-context.
     * @param key the name of the attribute containing the double array.
     * @param defaultValue the double array to return in case there is no such an annotation.
     * @param missingValue value to use to fill up positions with a missing value (e.g. '.').
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static int[] getAttributeAsIntArray(final Genotype genotype, final String key,
                                                     final Supplier<int[]> defaultValue, final int missingValue) {
        Utils.nonNull(genotype);
        return attributeValueToIntArray(genotype.getExtendedAttribute(key), key, defaultValue, missingValue);
    }

    private static double[] attributeValueToDoubleArray(final Object value, final String key, final Supplier<double[]> defaultResult, final double missingValue) {
        Utils.nonNull(key);
        final ToDoubleFunction<Object> doubleConverter = o -> {
            if (o == null) {
                return missingValue;
            } else {
                final String s = String.valueOf(o);
                if (s.equals(VCFConstants.MISSING_VALUE_v4)) {
                    return missingValue;
                } else {
                    try {
                        return Double.parseDouble(s);
                    } catch (final NumberFormatException ex) {
                        throw new GATKException(String.format("INFO annotation '%s' contains a non-double value '%s'", key, s), ex);
                    }
                }
            }
        };

        if (value == null) {
            return defaultResult.get();
        } else if (value.getClass().isArray()) {
            final double[] result = new double[Array.getLength(value)];
            for (int i = 0; i < result.length; i++) {
                result[i] = doubleConverter.applyAsDouble(String.valueOf(Array.get(value, i)));
            }
            return result;
        } else if (value.getClass().isAssignableFrom(Iterable.class)) {
            return StreamSupport.stream(((Iterable<?>)value).spliterator(), false)
                    .mapToDouble(doubleConverter).toArray();
        } else { // as a last resort with transform it into an String and try to parse an array out of it.
            return Stream.of(String.valueOf(value).trim().replaceAll("\\[|\\]", "")
                    .split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR))
                    .mapToDouble(doubleConverter).toArray();
        }
    }

    private static int[] attributeValueToIntArray(final Object value, final String key, final Supplier<int[]> defaultResult, final int missingValue) {
        Utils.nonNull(key);
        final ToIntFunction<Object> intConverter = o -> {
            if (o == null) {
                return missingValue;
            } else {
                final String s = String.valueOf(o).trim();
                if (s.equals(VCFConstants.MISSING_VALUE_v4)) {
                    return missingValue;
                } else {
                    try {
                        return Integer.parseInt(s);
                    } catch (final NumberFormatException ex) {
                        throw new GATKException(String.format("INFO annotation '%s' contains a non-int value '%s'", key, s), ex);
                    }
                }
            }
        };

        if (value == null) {
            return defaultResult.get();
        } else if (value.getClass().isArray()) {
            final int[] result = new int[Array.getLength(value)];
            for (int i = 0; i < result.length; i++) {
                result[i] = intConverter.applyAsInt(String.valueOf(Array.get(value, i)));
            }
            return result;
        } else if (value.getClass().isAssignableFrom(Iterable.class)) {
            return StreamSupport.stream(((Iterable<?>)value).spliterator(), false)
                    .mapToInt(intConverter).toArray();
        } else { // as a last resort with transform it into an String and try to parse an array out of it.
            return Stream.of(String.valueOf(value).trim().replaceAll("\\[|\\]", "")
                    .split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR))
                    .mapToInt(intConverter).toArray();
        }
    }

    /**
     * Returns an attribute as a string.
     * @param genotype the source genotype.
     * @param key the attribute key.
     * @param defaultValue default value to return in case the attribute does not have a value defined.
     * @return the string version of the attribute value.
     * @throws IllegalArgumentException if {@code genotype} or {@code key} is {@code null}.
     */
    public static String getAttributeAsString(final Genotype genotype, final String key, final String defaultValue) {
        Utils.nonNull(genotype);
        Utils.nonNull(key);
        final Object value = genotype.getExtendedAttribute(key);
        if (value == null) {
            return defaultValue;
        } else {
            return String.valueOf(value);
        }
    }

    /**
     * Returns an attribute as an integer.
     * @param genotype the source genotype.
     * @param key the attribute key.
     * @param defaultValue default value to return in case the attribute does not have a value defined.
     * @return the string version of the attribute value.
     * @throws IllegalArgumentException if {@code genotype} or {@code key} is {@code null}.
     * @throws NumberFormatException if the annotation exists, and its string transformation cannot be convert into an integer.
     */
    public static int getAttributeAsInt(final Genotype genotype, final String key, final int defaultValue) {
        Utils.nonNull(genotype);
        Utils.nonNull(key);
        final Object value = genotype.getExtendedAttribute(key);
        if (value == null) {
            return defaultValue;
        } else {
            try {
                return Integer.parseInt(String.valueOf(value));
            } catch (final NumberFormatException ex) {
                throw new NumberFormatException(String.format("attribute '%s' does not have a valid integer value: '%s'", key, String.valueOf(value)));
            }
        }
    }

    /**
     * Returns an attribute as a double.
     * @param genotype the source genotype.
     * @param key the attribute key.
     * @param defaultValue default value to return in case the attribute does not have a value defined.
     * @return the string version of the attribute value.
     * @throws IllegalArgumentException if {@code genotype} or {@code key} is {@code null}.
     * @throws NumberFormatException if the annotation exists, and its string transformation cannot be convert into a double.
     */
    public static double getAttributeAsDouble(final Genotype genotype, final String key, final double defaultValue) {
        Utils.nonNull(genotype);
        Utils.nonNull(key);
        final Object value = genotype.getExtendedAttribute(key);
        if (value == null) {
            return defaultValue;
        } else {
            try {
                return Double.parseDouble(String.valueOf(value));
            } catch (final NumberFormatException ex) {
                throw new NumberFormatException(String.format("attribute '%s' does not have a valid double value: '%s'", key, String.valueOf(value)));
            }
        }
    }

    /**
     * Instead of using the GQ value, it re-calculates the quality from the PL so that it does not need
     * to be bounded by an artificial maximum such as the standard GQ = 99.
     * @param genotype the PL source genotype.
     *
     * @return never {@code null}, but {@link Double#NaN} if there is no PLs.
     * @throws IllegalArgumentException if {@code genotype} is {@code null}
     */
    public static double calculateGenotypeQualityFromPLs(final Genotype genotype) {
        Utils.nonNull(genotype);
        if (!genotype.hasPL()) {
            return Double.NaN;
        } else {
            final int GQ = MathUtils.secondSmallestMinusSmallest(genotype.getPL(), -1);
            return GQ < 0 ? Double.NaN : (double) GQ;
        }
    }

    /**
     * Instead of using the GQ value, it re-calculates the quality from the PL so that it does not need
     * to be bounded by an artificial maximum such as the standard GQ = 99.
     * @param builder where to set the genotypes.
     * @param genotype the PL source genotype.
     *
     * @throws IllegalArgumentException if {@code genotype} is {@code null}
     */
    public static void setGenotypeQualityFromPLs(final GenotypeBuilder builder, final Genotype genotype) {
        final double gq = calculateGenotypeQualityFromPLs(genotype);
        if (Double.isNaN(gq)) {
            builder.noGQ();
        } else {
            builder.GQ((int) Math.floor(gq));
        }
    }

    /**
     * Gets an attribute value by transforming the string its string representation using an translation function.
     * @param g the source genotype.
     * @param key the attribute key.
     * @param translate the function to translate from an string to the return class of interest.
     * @param defaultValue the default value to return in case the attribute is not defined or is declared as missing.
     * @param <T> the return object type.
     * @return {@code null} or an instance of the class {@link T}
     */
    public static <T> T getAttributeAsObject(final Genotype g, final String key, final Function<String, T> translate,
                                                                              final T defaultValue) {
        Utils.nonNull(g);
        Utils.nonNull(key);
        final Object value = g.getExtendedAttribute(key);
        if (value == null || VCFConstants.MISSING_VALUE_v4.equals(value)) {
            return defaultValue;
        } else {
            return translate.apply(String.valueOf(value));
        }
    }

    /**
     * Get the pileup of reads covering a locus.  This is useful, for example, in VariantWalkers, which work on
     * ReadsContexts and not AlignmentContexts.
     */
    public static ReadPileup getPileup(final Locatable loc, final Iterable<GATKRead> reads) {
        final List<PileupElement> pile = StreamSupport.stream(reads.spliterator(), false)
                .filter(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK.and(ReadFilterLibrary.NOT_DUPLICATE))
                .map(AlignmentStateMachine::new)
                .map(asm -> {
                    while ( asm.stepForwardOnGenome() != null && asm.getGenomePosition() < loc.getStart()) { }
                    return asm.getGenomePosition() == loc.getStart() ? asm.makePileupElement() : null;
                }).filter(Objects::nonNull).collect(Collectors.toList());

        return new ReadPileup(loc, pile);
    }

    /** This is lifted directly from htsjdk with some minor modifications!  However, it is a private method there.
     *
     * This method cannot return {@link VariantContext.Type} MIXED
     *
     * Please see https://github.com/samtools/htsjdk/issues/999
     *
     * <p>Here are some cases that will not work properly, though this may not be an issue in practice:  </p>
     *  <ul>
     *      <li>"CGT" --> "GGA" this will be a MNP, but really it is two SNPs.</li>
     *      <li>Spanning deletions for alternate will show as {@link VariantContext.Type} NO_VARIATION</li>
     *      <li>Spanning deletions for reference will throw exception. </li>
     *      <li>Reference that is symbolic will throw an exception.</li>
     *  </ul>
     *
     * @param ref reference allele. Never {@code null}
     * @param allele alternate allele to compare. Never {@code null}
     * @return
     */
    public static VariantContext.Type typeOfVariant(final Allele ref, final Allele allele) {
        Utils.nonNull(ref);
        Utils.nonNull(allele);

        if ( ref.isSymbolic() )
            throw new IllegalStateException("Unexpected error: encountered a record with a symbolic reference allele");

        if ( allele.isSymbolic() )
            return VariantContext.Type.SYMBOLIC;

        if (allele.equals(Allele.SPAN_DEL)) {
            return VariantContext.Type.NO_VARIATION;
        }

        if ( ref.equals(Allele.SPAN_DEL) )
            throw new IllegalStateException("Unexpected error: encountered a record with a spanning deletion reference allele");

        if ( ref.length() == allele.length() ) {
            if (ref.basesMatch(allele)) {
                return VariantContext.Type.NO_VARIATION;
            } else if ( allele.length() == 1 )
                return VariantContext.Type.SNP;

            // If the two alleles are the same length and only differ by one base, then still a SNP.
            else if (IntStream.range(0, ref.length()).filter(i -> ref.getBases()[i] != allele.getBases()[i]).count() == 1) {
                return VariantContext.Type.SNP;
            } else
                return VariantContext.Type.MNP;
        }

        // Important note: previously we were checking that one allele is the prefix of the other.  However, that's not an
        // appropriate check as can be seen from the following example:
        // REF = CTTA and ALT = C,CT,CA
        // This should be assigned the INDEL type but was being marked as a MIXED type because of the prefix check.
        // In truth, it should be absolutely impossible to return a MIXED type from this method because it simply
        // performs a pairwise comparison of a single alternate allele against the reference allele (whereas the MIXED type
        // is reserved for cases of multiple alternate alleles of different types).  Therefore, if we've reached this point
        // in the code (so we're not a SNP, MNP, or symbolic allele), we absolutely must be an INDEL.

        return VariantContext.Type.INDEL;

        // old incorrect logic:
        // if (oneIsPrefixOfOther(ref, allele))
        //     return Type.INDEL;
        // else
        //     return Type.MIXED;
    }

    /**
     *  This method should only be run on variants that are known to be indels.  See {@link GATKProtectedVariantContextUtils::typeOfVariant}
     *
     *<p>Here are some cases that will not work properly, though this may not be an issue in practice:  </p>
     *  <ul>
     *      <li>"CT" --> "CATT" this is really just a simple AT insertion, but this will show up as complex.</li>
     *  </ul>
     * @param ref reference allele. Never {@code null}
     * @param allele alternate allele to compare. Never {@code null}
     * @return true if the indel is complex (for example, also includes a SNP), false if simple indel.  If the input alleles define a variant that is not
     *  an indel, then the behavior of this method is undefined (though will probably just return false).
     *
     */
    public static boolean isComplexIndel(final Allele ref, final Allele allele) {

        Utils.nonNull(ref);
        Utils.nonNull(allele);

        // Symbolic --> false
        if (ref.isSymbolic() || (ref.length() == 0)) {
            return false;
        }
        if (allele.isSymbolic() || (allele.length() == 0)) {
            return false;
        }

        // SNP, MNP, or no variation --> false
        if ( ref.length() == allele.length() ) {
            return false;
        }

        // obvious simple del or simple indel
        if ((allele.length() == 1) || (ref.length() == 1)) {
            return false;
        }

        // If the ref starts with the alt or vice versa, this is still simple.
        if (allele.length() > ref.length()) {
            final boolean isAltStartsWithRef = IntStream.range(0, ref.length()).allMatch(i -> ref.getBases()[i] == allele.getBases()[i]);
            return !isAltStartsWithRef;
        } else {
            final boolean isRefStartsWithAlt = IntStream.range(0, allele.length()).allMatch(i -> ref.getBases()[i] == allele.getBases()[i]);
            return !isRefStartsWithAlt;
        }
    }

    /**
     * Given a set of alleles (reference and alternate), choose the allele that is the best match for the given read (and offset)
     * TODO: This method cannot recognize equivalent alleles (See https://github.com/broadinstitute/gatk/issues/5061)
     * @param pileupElement read and offset.  Never {@code null}
     * @param referenceAllele Reference allele.  Never {@code null}
     * @param altAlleles List of candidate alternate alleles.  Never {@code null}
     * @param minBaseQualityCutoff minimum base quality for the bases that match the allele in order to be counted.
     *                             Must be positive or zero.  If you do not want any filtering, specify 0.
     * @return The allele (reference or from the altAlleles) that matches.  {@code null} if none are a match or the base qualities
     *      corresponding to the allele don't all exceed the minimum.
     */
    public static Allele chooseAlleleForRead(final PileupElement pileupElement, final Allele referenceAllele, final List<Allele> altAlleles, int minBaseQualityCutoff) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(referenceAllele);
        Utils.nonNull(altAlleles);
        ParamUtils.isPositiveOrZero(minBaseQualityCutoff, "Minimum base quality must be positive or zero.");

        final boolean isRef = referenceAllele.basesMatch(getBasesForAlleleInRead(pileupElement, referenceAllele))
                && !pileupElement.isBeforeDeletionStart() && !pileupElement.isBeforeInsertion();

        Allele pileupAllele = null;
        if (!isRef) {

            for (Allele altAllele : altAlleles) {
                final VariantContext.Type variantType = typeOfVariant(referenceAllele, altAllele);

                if (variantType == VariantContext.Type.INDEL) {
                    if (isIndelInThePileupElement(pileupElement, referenceAllele, altAllele)) {
                        pileupAllele = altAllele;
                    }

                } else if (variantType == VariantContext.Type.MNP || variantType == VariantContext.Type.SNP) {
                    if (doesReadContainAllele(pileupElement, altAllele) == Trilean.TRUE) {
                        pileupAllele = altAllele;
                    }
                }
            }
        } else {
            pileupAllele = referenceAllele;
        }

        if ((pileupAllele != null) && (getMinBaseQualityForAlleleInRead(pileupElement, pileupAllele) < minBaseQualityCutoff)) {
            pileupAllele = null;
        }

        return pileupAllele;
    }

    private static boolean isIndelInThePileupElement(final PileupElement pileupElement, final Allele referenceAllele, final Allele altAllele) {
        boolean isAltAlleleInThePileup = false;

        // Check insertion
        if (pileupElement.isBeforeInsertion()) {
            final String insertionBases = pileupElement.getBasesOfImmediatelyFollowingInsertion();
            // edge case: ignore a deletion immediately preceding an insertion as p.getBasesOfImmediatelyFollowingInsertion() returns null [EB]
            if ((insertionBases != null) && (Allele.extend(referenceAllele, insertionBases.getBytes()).basesMatch(altAllele))) {
                isAltAlleleInThePileup = true;
            }
        } else if (pileupElement.isBeforeDeletionStart()) {
            final int deletionLength = pileupElement.getLengthOfImmediatelyFollowingIndel();
            if ((referenceAllele.getBases().length - altAllele.getBases().length) == deletionLength) {
                isAltAlleleInThePileup = true;
            }
        }
        return isAltAlleleInThePileup;
    }

    /**
     * @param pileupElement pileup element representing the read.  Never {@code null}
     * @param allele allele to get overlapping bases in read.  Never {@code null}
     * @return array of the bytes that correspond to the allele in the pileup element.  Note that, if the read ends, this
     * list can be smaller than the length of the allele.
     */
    private static byte[] getBasesForAlleleInRead(final PileupElement pileupElement, final Allele allele) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(allele);
        return ArrayUtils.subarray(pileupElement.getRead().getBases(), pileupElement.getOffset(), pileupElement.getOffset() + allele.getBases().length);
    }

    /**
     * TODO: Test.  And make sure to test with reference alleles to make sure "*" is not included.
     * @param pileupElement pileup element representing the read.  Never {@code null}
     * @param allele query allele.  Never {@code null}
     * @return Whether the read contains the allele.  Note that unknown can occur as well.
     */
    public static Trilean doesReadContainAllele(final PileupElement pileupElement, final Allele allele) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(allele);

        final byte[] readBases = ArrayUtils.subarray(pileupElement.getRead().getBases(), pileupElement.getOffset(), pileupElement.getOffset() + allele.getBases().length);

        if (readBases.length < allele.getBases().length) {
            return Trilean.UNKNOWN;
        }

        if (allele.basesMatch(readBases)) {
            return Trilean.TRUE;
        } else {
            return Trilean.FALSE;
        }
    }

    /**
     * Find the minimum base quality for all bases in a read that correspond to a given allele.
     *
     * @param pileupElement pileup element representing the read.  Never {@code null}
     * @param allele query allele.  Never {@code null}
     * @return lowest base quality seen in the corresponding bases of the read
     */
    private static int getMinBaseQualityForAlleleInRead(final PileupElement pileupElement, final Allele allele) {
        Utils.nonNull(pileupElement);
        Utils.nonNull(allele);
        final byte[] alleleBases = allele.getBases();
        final byte[] pileupBaseQualities = ArrayUtils.subarray(pileupElement.getRead().getBaseQualities(), pileupElement.getOffset(), pileupElement.getOffset() + alleleBases.length);
        final OptionalInt minQuality = IntStream.range(0, pileupBaseQualities.length).map(i -> Byte.toUnsignedInt(pileupBaseQualities[i])).min();
        return minQuality.orElse(-1);
    }
}
