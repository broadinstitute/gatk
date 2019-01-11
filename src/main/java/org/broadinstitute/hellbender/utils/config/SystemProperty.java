package org.broadinstitute.hellbender.utils.config;

import java.lang.annotation.*;

/**
 * An annotation to denote Configuration options that should be injected into the Java System Properties.
 * The presence of this annotation on a {@link org.aeonbits.owner.Config} object marks it for addition into the
 * Java System Properties in {@link ConfigFactory#injectSystemPropertiesFromConfig}.
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.METHOD)
@Documented
public @interface SystemProperty {

}
