package org.broadinstitute.hellbender.cmdline;

import joptsimple.ValueConversionException;
import joptsimple.ValueConverter;

/**
 * converts values case insensitively matching T, True, F, or False to true or false
 * throws {@link ValueConversionException} otherwise
 */
public final class StrictBooleanConverter implements ValueConverter<String> {
    public String convert( String value ) {
        if ( value.equalsIgnoreCase("true") || value.equalsIgnoreCase("t")) {
            return "true";
        } else if (value.equalsIgnoreCase("false") || value.equalsIgnoreCase("f")) {
            return "false";
        } else {
            throw new ValueConversionException(value + " does not match one of T|True|F|False");
        }
    }
    public final Class<? extends String> valueType() {
        return String.class;
    }

    public String valuePattern() {
        return "[T|True|F|False]";
    }
}
