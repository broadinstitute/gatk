package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;

import java.util.LinkedHashSet;
import java.util.Map;

/**
 * Abstract class representing a {@link Funcotator} annotation.
 * Created by jonn on 8/30/17.
 */
public interface Funcotation {

    /**
     * @return The alternate {@link Allele} that is associated with this {@link Funcotation}.
     */
    Allele getAltAllele();

    /**
     * Override the given field with a given value for when it comes time to serialize and write this {@link Funcotation}.
     * If the given {@code field} is not contained in this {@link Funcotation} then a {@link org.broadinstitute.hellbender.exceptions.UserException} will be thrown.
     * @param fieldName A {@link String} comprising the name of the field to override.
     * @param overrideValue The {@link String} value for the field to override.
     */
    void setFieldSerializationOverrideValue( final String fieldName, final String overrideValue );

    /**
     * Converts this {@link Funcotation} to a string suitable for insertion into a VCF file.
     * @return a {@link String} representing this {@link Funcotation} suitable for insertion into a VCF file.
     */
    String serializeToVcfString();

    /**
     * @return The name of the data source behind the {@link DataSourceFuncotationFactory} used to create this {@link Funcotation}.
     */
    String getDataSourceName();

    /**
     * Get the names of the fields in this {@link Funcotation}.
     * @return The ordered set of fields in this {@link Funcotation} as a {@link LinkedHashSet} of {@link String}s.
     */
    LinkedHashSet<String> getFieldNames();

    /**
     * Get the value of a field in this {@link Funcotation}.
     * @return The {@link String} value of a field in this {@link Funcotation}.
     * @throws {@link org.broadinstitute.hellbender.exceptions.GATKException} if the given {@code fieldName} is not in this {@link Funcotation}.
     */
    String getField(final String fieldName);

    /**
     * Get the value of a field in this {@link Funcotation}.
     * @return The {@link String} value of a field in this {@link Funcotation}.  If the field name is not present,
     * return the default value.
     */
    default String getFieldOrDefault(final String fieldName, final String defaultValue){
        if (hasField(fieldName)) {
            return getField(fieldName);
        } else {
            return defaultValue;
        }
    }

    /**
     * Override fields with values as specified by the input map (for when it comes time to serialize and write this {@link Funcotation}).
     * If the given overrides map is null, will not override any field.
     * If the given overrides map is not null and if the given {@code field} is not contained in this {@link Funcotation} then a {@link org.broadinstitute.hellbender.exceptions.UserException} will be thrown.
     * @param fieldSerializationOverrides A {@link Map} containing fields to override in this {@link Funcotation}.
     */
    default void setFieldSerializationOverrideValues(final Map<String,String> fieldSerializationOverrides) {
        if (fieldSerializationOverrides != null) {
            for ( final String field : fieldSerializationOverrides.keySet() ) {
                setFieldSerializationOverrideValue(field, fieldSerializationOverrides.get(field));
            }
        }
    }

    /**
     * TODO: This interface should have nothing specific to a VCF.  That should be the job of the VCFOutputRenderer to sanitize any strings.  https://github.com/broadinstitute/gatk/issues/4797
     * Converts this {@link Funcotation} to a string suitable for insertion into a VCF file.
     * {@code manualAnnotationString} should be written first, followed by the inherent annotations in this {@link Funcotation}.
     * @param manualAnnotationString A {@link String} of manually-provided annotations to add to this {@link Funcotation}.
     * @return a {@link String} representing this {@link Funcotation} suitable for insertion into a VCF file.
     */
    default String serializeToVcfString(final String manualAnnotationString) {
        return (manualAnnotationString == null ? "" : manualAnnotationString) + serializeToVcfString();
    }

    /**
     * @return Return whether the field exists in this {@link Funcotation}.
     */
    boolean hasField(final String fieldName);

    /**
     * @return Metadata for this {@link Funcotation}.  Never {@code null}.  All fields in {@see getFieldNames} should be
     * represented.
     */
    FuncotationMetadata getMetadata();
}
