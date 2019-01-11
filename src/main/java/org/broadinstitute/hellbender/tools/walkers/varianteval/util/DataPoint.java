package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface DataPoint {
    String description() default ""; // the description, optional
    String format() default "";
}
