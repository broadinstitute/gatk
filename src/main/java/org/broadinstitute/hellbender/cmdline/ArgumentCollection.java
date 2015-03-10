package org.broadinstitute.hellbender.cmdline;

import java.lang.annotation.*;

/**
 * Used to annotate a field in a CommandLineProgram that holds a instance containing @Option-annotated
 * fields.  To set a value for a nested option on the command line, use <member-name>.<option>=value.
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
@Documented
@Inherited
public @interface ArgumentCollection {
    /** Text that appears for this group of options in text describing usage of the command line program. */
    String doc() default "";
}
