package org.broadinstitute.hellbender.cmdline;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Annotating a command line tool with this interface targets it for WDL generation. Any tool on which this annotation
 * is used must have one or more output files of type {@link org.broadinstitute.hellbender.engine.GATKOutputPath}.
 */
@Target(ElementType.TYPE)
@Retention(RetentionPolicy.RUNTIME)
public @interface RuntimeProperties {

    /**
     * @return a WDL-compatible string specifying the runtime memory requirements for this tool
     */
    String memory() default "";
}
