package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;

/**
 * An abstract class to allow for writing output for the Funcotator.
 * Used to output Funcotations to a location of the user's choice.
 * For example, writing out to a VCF file (e.g. {@link org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer}).
 * It is also possible to write to a stream, or other connection (e.g. a UDP port)
 * at the user's discretion.
 * Created by jonn on 8/30/17.
 */
public interface OutputRenderer extends AutoCloseable {

    /**
     * Open the {@link OutputRenderer} for writing.
     */
    void open();

    /**
     * Close the {@link OutputRenderer}.
     */
    void close();

    /**
     * Write the given {@code variant} and {@code funcotations} to the output file.
     * @param variant {@link VariantContext} to write to the file.
     * @param funcotations {@link List} of {@link Funcotation} to add to the given {@code variant} on output.
     */
    void write(final VariantContext variant, final List<Funcotation> funcotations);
}
