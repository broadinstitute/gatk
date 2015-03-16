package org.broadinstitute.hellbender.cmdline;

import java.io.Serializable;

/**
 * marker interface for classes that are intended to be used with @ArgumentCollection
 * Those are parsed by the CommandLineParser class.
 */
public interface ArgumentCollectionDefinition extends Serializable{
}
