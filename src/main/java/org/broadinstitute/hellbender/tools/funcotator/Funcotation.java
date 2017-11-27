package org.broadinstitute.hellbender.tools.funcotator;

import java.util.Map;

/**
 * Abstract class representing a {@link Funcotator} annotation.
 * Created by jonn on 8/30/17.
 */
public abstract class Funcotation {

    /**
     * Override the given field with a given value for when it comes time to serialize and write this {@link Funcotation}.
     * If the given {@code field} is not contained in this {@link Funcotation} then a {@link org.broadinstitute.hellbender.exceptions.UserException} will be thrown.
     * @param field A {@link String} comprising the name of the field to override.
     * @param overrideValue The {@link String} value for the field to override.
     */
    public abstract void setFieldSerializationOverrideValue( final String field, final String overrideValue );

    /**
     * Override fields with values as specified by the input map (for when it comes time to serialize and write this {@link Funcotation}).
     * If the given {@code field} is not contained in this {@link Funcotation} then a {@link org.broadinstitute.hellbender.exceptions.UserException} will be thrown.
     * @param fieldSerializationOverrides A {@link Map} containing fields to override in this {@link Funcotation}.
     */
    public void setFieldSerializationOverrideValues(final Map<String,String> fieldSerializationOverrides) {
        for ( final String field : fieldSerializationOverrides.keySet() ) {
            setFieldSerializationOverrideValue(field, fieldSerializationOverrides.get(field));
        }
    }

    /**
     * Converts this {@link Funcotation} to a string suitable for insertion into a VCF file.
     * @return a {@link String} representing this {@link Funcotation} suitable for insertion into a VCF file.
     */
    public abstract String serializeToVcfString();

    /**
     * Converts this {@link Funcotation} to a string suitable for insertion into a VCF file.
     * {@code manualAnnotationString} should be written first, followed by the inherent annotations in this {@link Funcotation}.
     * @param manualAnnotationString A {@link String} of manually-provided annotations to add to this {@link Funcotation}.
     * @return a {@link String} representing this {@link Funcotation} suitable for insertion into a VCF file.
     */
    public abstract String serializeToVcfString(final String manualAnnotationString);
}
