package org.broadinstitute.hellbender.cmdline;

import java.lang.annotation.*;

/**
 * Annotates a command line program with various properties, such as usage (short and long),
 * as well as to which program group it belongs.
 *
 * TODO: enforced that any CommandLineProgram has this property defined (use an annotation processor?).
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
@Inherited
public @interface CommandLineProgramProperties {
    /**
     * @return a summary of what the program does
     */
    String summary();

    /**
     * @return a very short summary for the main menu list of all programs
     */
    String oneLineSummary();

    /**
     * @return an example command line for this program
     */
    String usageExample() default "The author of this program hasn't included any example usage, please complain to them.";
    Class<? extends CommandLineProgramGroup> programGroup();
    boolean omitFromCommandLine() default false;
}
