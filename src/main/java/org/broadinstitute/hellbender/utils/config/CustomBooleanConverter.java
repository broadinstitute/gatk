package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Converter;

import java.lang.reflect.Method;

/**
 * Converts a given string into a Boolean after trimming whitespace from that string.
 * Allows a boolean value with a trailing space to be correctly parsed.
 * This exists to show how to create custom converter classes for types when reading from a config file with {@link org.aeonbits.owner}.
 */
final public class CustomBooleanConverter implements Converter<Boolean> {

    @Override
    public Boolean convert(final Method method, final String input) {
        return Boolean.parseBoolean( input.trim() );
    }
}
