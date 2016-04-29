package org.broadinstitute.hellbender.cmdline;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Used to annotate which field of a CommandLineProgram should store parameters given at the 
 * command line which are not options. Fields with this annotation must be a Collection
 * (and probably should be a List if order is important).
 * If a command line call looks like "cmd option=foo x=y bar baz" the values "bar" and "baz"
 * would be added to the collection with this annotation. The java type of the arguments
 * will be inferred from the generic type of the collection. The type must be an enum or
 * have a constructor with a single String parameter.
 *
 * @author Alec Wysoker
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
@Documented
public @interface PositionalArguments {
    /** The minimum number of arguments required. */
    int minElements() default 0;
    
    /** The maximum number of arguments allowed. */
    int maxElements() default Integer.MAX_VALUE;

    /**
     * Documentation for the command-line argument.  Should appear when the
     * --help argument is specified.
     * @return Doc string associated with this command-line argument.
     */
    String doc() default "Undocumented option";
}
