package org.broadinstitute.hellbender.cmdline;

import java.lang.annotation.*;

/**
 * Indicates that a walker or walker argument should not be presented in the help system.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.FIELD})
public @interface Hidden {
}
