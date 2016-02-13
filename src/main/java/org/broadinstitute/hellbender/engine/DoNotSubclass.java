package org.broadinstitute.hellbender.engine;

import java.lang.annotation.*;

/**
 * Classes annotated with this annotation are NOT intended or designed to be extended and should be treated as final.
 * The reason to make them non-final is to enable testing using mock objects.
 */
@Documented
@Retention(RetentionPolicy.CLASS)
@Target(ElementType.TYPE)
public @interface DoNotSubclass {
}
