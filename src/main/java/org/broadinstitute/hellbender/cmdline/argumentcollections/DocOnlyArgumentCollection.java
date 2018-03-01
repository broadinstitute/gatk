package org.broadinstitute.hellbender.cmdline.argumentcollections;

/**
 * Marker interface for documentation-only arguments.

 * <p>Some arguments might be exposed on tools only for documentation, even if they are parsed in Main. Possible use
 * cases include:
 *
 * <ul>
 *     <li>Wrapper script arguments that might appear in the command line but are unused.</li>
 *      <li>Arguments that should be parsed on Main to initialize custom configurations (e.g., defaults loaded
 *          as static initializers.</li>
 * </ul>
 *
 * <p>This interface allows to override documentation-only arguments to allow downstream projects to expose/hide
 * specific documentation arguments of GATK or include custom ones.
 *
 * @author Daniel Gomez-Sanchez (magicDGS)
 * @see GATKDocOnlyArgumentCollection for an example implementation.
 */
public interface DocOnlyArgumentCollection {

}
