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

    // The arguments in this collection are dependent on (and only enabled with the presence of)
    // another argument. This may be used alone or together with dependsOnValue.
    String dependsOnArgument() default "";

    // The arguments in this collection are dependent on (and only enabled with the presence of)
    // another argument with this specific value.
    String dependsOnValue() default "";
}
