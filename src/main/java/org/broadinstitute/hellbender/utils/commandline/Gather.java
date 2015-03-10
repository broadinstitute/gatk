package org.broadinstitute.hellbender.utils.commandline;

import java.lang.annotation.*;

/**
 * Specifies the class type to gather an @Output
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.FIELD})
@SuppressWarnings("rawtypes")
public @interface Gather {
    Class value() default Gather.class;
    String className() default "";
    boolean enabled() default true;
}
