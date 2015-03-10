package org.broadinstitute.hellbender.utils.commandline;

import java.lang.annotation.*;

/**
 * Indicates that a walker or walker argument should not be presented in the help system.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.TYPE, ElementType.FIELD})
public @interface HiddenOption {
}
