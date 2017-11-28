package org.broadinstitute.hellbender.tools.funcotator.dataSources.XSV;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A {@link Funcotation} to hold data from XSV files.
 * Created by jonn on 11/28/17.
 */
public class XSVFuncotation extends Funcotation {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    /**
     * Names of the fields in this {@link XSVFuncotation}.
     */
    private LinkedHashMap<String, String> fieldMap;

    //==================================================================================================================
    // Constructors:

    public XSVFuncotation( final List<String> fieldNames, final List<String> fieldValues ) {
        if ( fieldNames.size() != fieldValues.size() ) {
            throw new UserException.BadInput("Field names and Field values are of different lengths!  This must not be!");
        }

        fieldMap = new LinkedHashMap<>(fieldNames.size());
        for ( int i = 0; i < fieldNames.size() ; ++i ) {
            fieldMap.put(fieldNames.get(i), fieldValues.get(i));
        }
    }

    public XSVFuncotation( final LinkedHashMap<String, String> fieldMap ) {
        this.fieldMap = fieldMap;
    }

    public XSVFuncotation( final XSVFuncotation that ) {
        this.fieldMap = new LinkedHashMap<>( that.fieldMap );
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void setFieldSerializationOverrideValue(final String fieldName, final String overrideValue) {
        if ( !fieldMap.containsKey(fieldName) ) {
            throw new GATKException("Attempted to override a field that is not contained in this XSVFuncotation: "
                    + fieldName + " is not one of [" + String.join(",", fieldMap.keySet()) + "]");
        }
        fieldMap.put(fieldName, overrideValue);
    }

    @Override
    public String serializeToVcfString() {
        return fieldMap.values().stream()
                .map(f -> (f == null ? "" : f))
                .collect(Collectors.joining(VcfOutputRenderer.FIELD_DELIMITER));
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:
}
