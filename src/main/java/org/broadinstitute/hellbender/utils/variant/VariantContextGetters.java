package org.broadinstitute.hellbender.utils.variant;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class VariantContextGetters {

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
     * Composes the int array from an INFO annotation.
     *
     * @param variantContext the target variant-context.
     * @param key the name of the attribute containing the double array.
     * @param defaultValue the double array to return in case there is no such an annotation.
     * @param missingValue value to use to fill up positions with a missing value (e.g. '.').
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static int[] getAttributeAsIntArray(final VariantContext variantContext, final String key,
                                               final Supplier<int[]> defaultValue, final int missingValue) {
        Utils.nonNull(variantContext);
        return attributeValueToIntArray(variantContext.getAttribute(key), key, defaultValue, missingValue);
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
     * Composes the double array from a genotype annotation. Provides default and missing values.
     *
     * @param variantContext the target variant-context.
     * @param attribute the name of the attribute containing the double array.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static double[] getAttributeAsDoubleArray(final VariantContext variantContext, final String attribute) {
        return getAttributeAsDoubleArray(variantContext, attribute, () -> null, -1);
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
     * Get Long attribute from a variant context.
     *
     * @param variantContext the target variant-context.
     * @param attribute the name of the attribute containing the Long value.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static Long getAttributeAsLong(final VariantContext variantContext, final String attribute, final Long defaultValue) {
        Utils.nonNull(variantContext);
        Utils.nonNull(attribute);
        Object x = variantContext.getAttribute(attribute);
        if ( x == null || x.equals(VCFConstants.MISSING_VALUE_v4) ) return defaultValue;
        if ( x instanceof Number ) return ((Number) x).longValue();
        return Long.valueOf((String)x); // throws an exception if this isn't a string
    }

    /**
     * Composes the Long List from a variant context.
     *
     * @param variantContext the target variant-context.
     * @param attribute the name of the attribute containing the Long list.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static List<Long> getAttributeAsLongList(final VariantContext variantContext, final String attribute, final Long defaultValue) {
        Utils.nonNull(variantContext);
        Utils.nonNull(attribute);
        return variantContext.getAttributeAsList(attribute).stream().map(
                x -> {
                    if (x == null || x.equals(VCFConstants.MISSING_VALUE_v4)) {
                        return defaultValue;
                    } else if (x instanceof Number) {
                        return ((Number) x).longValue();
                    } else {
                        return Long.valueOf((String)x); // throws an exception if this isn't a string
                    }
                }
        ).collect(Collectors.toList());
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

    //copied from htsjdk.variant.variantcontext.CommonInfo.getAttributeAsList for simplicity
    //maybe we should expose this as a static method in htsjdk?
    @SuppressWarnings("unchecked")
    public static List<Object> attributeToList(final Object attribute){
        if ( attribute == null ) return Collections.emptyList();
        if ( attribute instanceof List) return (List<Object>)attribute;
        if ( attribute.getClass().isArray() ) {
            if (attribute instanceof int[]) {
                return Arrays.stream((int[])attribute).boxed().collect(Collectors.toList());
            } else if (attribute instanceof double[]) {
                return Arrays.stream((double[])attribute).boxed().collect(Collectors.toList());
            }
            return Arrays.asList((Object[])attribute);
        }
        if (attribute instanceof String) {
            return new ArrayList<>(Arrays.asList((Object[])((String)attribute).split(",")));
        }
        return Collections.singletonList(attribute);
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
}
