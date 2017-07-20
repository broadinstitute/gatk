package org.broadinstitute.hellbender.utils;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
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
            final int GQ = GATKProtectedMathUtils.secondSmallestMinusSmallest(genotype.getPL(), -1);
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
}
