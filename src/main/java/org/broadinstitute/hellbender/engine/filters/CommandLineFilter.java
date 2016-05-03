package org.broadinstitute.hellbender.engine.filters;

import java.io.Serializable;

/**
 * Interface implemented by command-line accessible read filters for argument validation.
 */
public interface CommandLineFilter extends Serializable {
    String validate();
}
