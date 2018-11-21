package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Molten for @Analysis modules.
 *
 * If you are flagged as a molten analysis, then there must be one and
 * only one annotation in that evaluation module: @Molten which
 * must have time Map<Object, Object>.  This data set will then
 * be represented in the VE output as:
 *
 * variable value
 * key1     value1
 * key2     value1
 * ...
 * keyN     valueN
 *
 * in the output table.  The names of these two fields can be override via annotation values.
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface Molten {
    String description() default ""; // the description, optional

    /**
     * The name to use for the molten variable field in the output table.
     * @return
     */
    String variableName() default "variable";
    String variableFormat() default "";

    /**
     * The name to use for the molten value field in the output table.
     * @return
     */
    String valueName() default "value";
    String valueFormat() default "";
}
