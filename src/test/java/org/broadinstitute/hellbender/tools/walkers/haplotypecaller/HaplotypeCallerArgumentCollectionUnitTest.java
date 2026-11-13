package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.map.DefaultedMap;
import org.apache.commons.lang3.reflect.FieldUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.util.*;
import java.util.stream.Collectors;

public class HaplotypeCallerArgumentCollectionUnitTest extends GATKBaseTest {

    // copyWithNewPloidy() makes the arguments for the genotyping engines of --ploidy-regions. Every argument, including
    // those of superclasses and nested argument collections, must be copied. The copy must share no argument collection
    // and no mutable argument value (List, Map, GATKPath or FeatureInput) with the original.
    @Test
    public void testArgCollectionCopiesAreSynchronized() throws IllegalAccessException {
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        final List<String> unhandled = new ArrayList<>();
        final int numChanged = setArgumentsToNonDefaultValues(hcArgs, "hcArgs", unhandled);
        // the per-sample contamination map (a DefaultedMap) is not an argument itself; it is set through its setter instead
        Assert.assertTrue(unhandled.remove("hcArgs.standardArgs.sampleContamination (DefaultedMap)"));
        hcArgs.standardArgs.setSampleContamination(Collections.singletonMap("sample", 0.1));
        Assert.assertTrue(unhandled.isEmpty(), "fields of types the test cannot change: " + unhandled);
        Assert.assertTrue(numChanged > 100, "only " + numChanged + " fields were changed");

        final int newPloidy = hcArgs.standardArgs.genotypeArgs.samplePloidy + 1;
        final HaplotypeCallerArgumentCollection hcArgsCopy = hcArgs.copyWithNewPloidy(newPloidy);
        Assert.assertEquals(hcArgsCopy.standardArgs.genotypeArgs.samplePloidy, newPloidy);
        Assert.assertEquals(hcArgs.standardArgs.genotypeArgs.samplePloidy, newPloidy - 1, "the original's ploidy changed");

        // apart from the ploidy, every field must be equal, and nothing mutable may be shared
        hcArgsCopy.standardArgs.genotypeArgs.samplePloidy = hcArgs.standardArgs.genotypeArgs.samplePloidy;
        final List<String> differences = new ArrayList<>();
        compareArguments(hcArgs, hcArgsCopy, "hcArgs", differences);
        Assert.assertTrue(differences.isEmpty(), "the copy differs from the original: " + differences);
    }

    /** All non-static fields of the object's class and its superclasses */
    private static List<Field> getInstanceFields(final Object object) {
        return FieldUtils.getAllFieldsList(object.getClass()).stream()
                .filter(field -> !Modifier.isStatic(field.getModifiers()))
                .peek(field -> field.setAccessible(true))
                .collect(Collectors.toList());
    }

    /** True for the GATK argument collections (and other GATK objects) nested in the arguments, which are compared field by field */
    private static boolean isNestedArgumentObject(final Class<?> type) {
        return type.getName().startsWith("org.broadinstitute") && !type.isEnum()
                && !GATKPath.class.isAssignableFrom(type);
    }

    /** True for argument values that can be changed in place, so the copy must not share them with the original */
    private static boolean isMutableValue(final Object value) {
        return value instanceof Collection || value instanceof Map || value instanceof GATKPath;
    }

    /**
     * Sets every field reachable from the object to a value other than its default, and returns how many fields were set.
     * Primitives, boxed types, Strings and enums are changed, Files, GATKPaths and FeatureInputs are replaced, and Lists get
     * one more element. Nested GATK objects are changed recursively. The paths of fields of any other type are added to unhandled.
     */
    private static int setArgumentsToNonDefaultValues(final Object object, final String path, final List<String> unhandled) throws IllegalAccessException {
        int numChanged = 0;
        for (final Field field : getInstanceFields(object)) {
            final String fieldPath = path + "." + field.getName();
            final Object value = field.get(object);
            if (value != null && isNestedArgumentObject(field.getType())) {
                numChanged += setArgumentsToNonDefaultValues(value, fieldPath, unhandled);
                continue;
            }
            final Object newValue = List.class.isAssignableFrom(field.getType()) ?
                    getListWithOneMoreElement(field, (List<?>) value, fieldPath) :
                    getNonDefaultValue(field.getType(), value, fieldPath);
            if (newValue == null) {
                unhandled.add(fieldPath + " (" + field.getType().getSimpleName() + ")");
            } else {
                field.set(object, newValue);
                numChanged++;
            }
        }
        return numChanged;
    }

    /** A copy of the list field's value with one more element, or null if the type of the elements is not handled */
    private static List<Object> getListWithOneMoreElement(final Field field, final List<?> value, final String path) {
        if (!(field.getGenericType() instanceof ParameterizedType)) {
            return null;
        }
        final Type elementType = ((ParameterizedType) field.getGenericType()).getActualTypeArguments()[0];
        final Type rawElementType = elementType instanceof ParameterizedType ? ((ParameterizedType) elementType).getRawType() : elementType;
        if (!(rawElementType instanceof Class)) {
            // e.g. a wildcard or type variable: reported as unhandled by the caller
            return null;
        }
        final Object newElement = getNonDefaultValue((Class<?>) rawElementType, null, path);
        if (newElement == null) {
            return null;
        }
        final List<Object> newList = value == null ? new ArrayList<>() : new ArrayList<>(value);
        newList.add(newElement);
        return newList;
    }

    /** A value of the type that differs from the given value, or null if the type is not handled */
    private static Object getNonDefaultValue(final Class<?> type, final Object value, final String path) {
        if (type == boolean.class || type == Boolean.class) {
            return value == null || !(Boolean) value;
        } else if (type == byte.class || type == Byte.class) {
            return (byte) (value == null ? 1 : (Byte) value + 1);
        } else if (type == int.class || type == Integer.class) {
            return value == null ? 1 : (Integer) value + 1;
        } else if (type == float.class || type == Float.class) {
            return value == null || !Float.isFinite((Float) value) ? 1.5f : (Float) value + 1.5f;
        } else if (type == double.class || type == Double.class) {
            return value == null || !Double.isFinite((Double) value) ? 1.5 : (Double) value + 1.5;
        } else if (type == String.class) {
            return "changed-" + path;
        } else if (type.isEnum()) {
            final Object[] constants = type.getEnumConstants();
            return constants[value == null ? 0 : (((Enum<?>) value).ordinal() + 1) % constants.length];
        } else if (type == File.class) {
            return new File(path + ".txt");
        } else if (type == GATKPath.class) {
            return new GATKPath(path + ".txt");
        } else if (type == FeatureInput.class) {
            return new FeatureInput<VariantContext>(new GATKPath(path + ".vcf"));
        }
        return null;
    }

    /**
     * Adds to differences the paths of fields that differ between the objects, and of nested GATK objects and mutable
     * values that the two share. DefaultedMaps must also have the same default value, which their equals() ignores.
     */
    private static void compareArguments(final Object original, final Object copy, final String path, final List<String> differences) throws IllegalAccessException {
        for (final Field field : getInstanceFields(original)) {
            final String fieldPath = path + "." + field.getName();
            final Object originalValue = field.get(original);
            final Object copyValue = field.get(copy);
            if (originalValue != null && originalValue == copyValue && (isNestedArgumentObject(field.getType()) || isMutableValue(originalValue))) {
                differences.add(fieldPath + ": the copy shares the original's instance");
            }
            if (originalValue != null && copyValue != null && isNestedArgumentObject(field.getType())) {
                compareArguments(originalValue, copyValue, fieldPath, differences);
            } else if (!Objects.equals(originalValue, copyValue)) {
                differences.add(fieldPath + ": " + originalValue + " in the original but " + copyValue + " in the copy");
            } else if (originalValue instanceof DefaultedMap && copyValue instanceof DefaultedMap) {
                final String absentKey = "key-not-in-" + fieldPath;
                final Object originalDefault = ((DefaultedMap<?, ?>) originalValue).get(absentKey);
                final Object copyDefault = ((DefaultedMap<?, ?>) copyValue).get(absentKey);
                if (!Objects.equals(originalDefault, copyDefault)) {
                    differences.add(fieldPath + ": default value " + originalDefault + " in the original but " + copyDefault + " in the copy");
                }
            }
        }
    }
}
