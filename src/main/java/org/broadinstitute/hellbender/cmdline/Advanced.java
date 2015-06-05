package org.broadinstitute.hellbender.cmdline;

import java.lang.annotation.*;

/**
 * Indicates that a walker argument should is considered an advanced option.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.FIELD})
public @interface Advanced {
}
