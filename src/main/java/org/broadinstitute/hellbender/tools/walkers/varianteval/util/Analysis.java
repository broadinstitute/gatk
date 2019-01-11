package org.broadinstitute.hellbender.tools.walkers.varianteval.util;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface Analysis {
    String name() default ""; // its description, required
    String description(); // its description, required
    boolean molten() default false; // if true we'll look for a @Molten map
}
