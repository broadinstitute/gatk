package org.broadinstitute.hellbender.tools.funcotator;

/**
 * Abstract class representing a {@link Funcotator} annotation.
 * Created by jonn on 8/30/17.
 */
public abstract class Funcotation {
    public abstract String serializeToVcfString();
}
