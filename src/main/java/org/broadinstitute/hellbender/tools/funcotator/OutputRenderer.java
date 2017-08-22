package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import java.util.List;

/**
 * An abstract class to allow for writing output for the Funcotator.
 * Created by jonn on 8/30/17.
 */
public abstract class OutputRenderer {

    /**
     * Open the {@link OutputRenderer} for writing.
     */
    public abstract void open();

    /**
     * Close the {@link OutputRenderer}.
     */
    public abstract void close();

    /**
     * Write the given {@code variant} and {@code funcotations} to the output file.
     * @param variant {@link VariantContext} to write to the file.
     * @param funcotations {@link List} of {@link Funcotation} to add to the given {@code variant} on output.
     */
    public abstract void write(final VariantContext variant, final List<Funcotation> funcotations);
}
