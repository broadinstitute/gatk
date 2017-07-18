package org.broadinstitute.hellbender.utils.config;

import org.aeonbits.owner.Converter;

import java.lang.reflect.Method;

/**
 * Created by jonn on 7/21/17.
 */
public class CustomBooleanConverter implements Converter<Boolean> {

    @Override
    public Boolean convert(Method method, String input) {
        return Boolean.getBoolean( input.trim() );
    }
}
