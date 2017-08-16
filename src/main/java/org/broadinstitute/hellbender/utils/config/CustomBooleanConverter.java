package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Converter;

import java.lang.reflect.Method;

/**
 * Converts a given string into a Boolean after trimming whitespace from that string.
 */
final public class CustomBooleanConverter implements Converter<Boolean> {

    @Override
    public Boolean convert(Method method, String input) {
        return Boolean.parseBoolean( input.trim() );
    }
}
