package org.broadinstitute.hellbender.cmdline;

import org.broadinstitute.hellbender.cmdline.programgroups.MiscProgramGroup;


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
    String usage();
    String usageShort();
    Class<? extends CommandLineProgramGroup> programGroup() default MiscProgramGroup.class;
    boolean omitFromCommandLine() default false;
}
