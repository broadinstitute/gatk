package org.broadinstitute.hellbender.utils.commandline;

import java.lang.annotation.*;

/**
 * Indicates that a walker argument should is considered an advanced option.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.TYPE,ElementType.FIELD})
public @interface AdvancedOption {
}
