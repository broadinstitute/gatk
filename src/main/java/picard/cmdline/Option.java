/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.cmdline;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Used to annotate which fields of a CommandLineProgram are options given at the command line.
 * If a command line call looks like "cmd option=foo x=y bar baz" the CommandLineProgram
 * would have annotations on fields to handle the values of option and x. All options
 * must be in the form name=value on the command line. The java type of the option
 * will be inferred from the type of the field or from the generic type of the collection
 * if this option is allowed more than once. The type must be an enum or
 * have a constructor with a single String parameter.
 *
 * @author Alec Wysoker
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.FIELD)
@Documented
public @interface Option {
    /** The name of the option as it would appear on the command line. */
    String shortName() default "";

    /** Text that appears for this option in text describing usage of the command line program. */
    String doc() default "";

    /**
     * If set to false, an exception will be thrown if the option is not specified.
     * If 2 options are mutually exclusive and both have optional=false it will be
     * interpreted as one or the other is required and an exception will only be thrown if
     * neither are specified.
     */
    boolean optional() default false;

    /**
     * Array of option names that cannot be used in conjunction with this one.
     * If 2 options are mutually exclusive and both have optional=false it will be
     * interpreted as one OR the other is required and an exception will only be thrown if
     * neither are specified.
     */
    String[] mutex() default {};

    /** The minimum number of times that this option is required. */
    int minElements() default 0;

    /** The maximum number of times this option is allowed. */
    int maxElements() default Integer.MAX_VALUE;

    /**
     * Is this an Option common to all command line programs.  If it is then it will only
     * be displayed in usage info when H or STDHELP is used to display usage.
     */
    boolean common() default false;

    /**
     * This boolean determines if this annotation overrides a parent annotation. If that is the case then
     * the options of the parent annotation are overridden with this annotation.
     */
    boolean overridable() default false;
}
