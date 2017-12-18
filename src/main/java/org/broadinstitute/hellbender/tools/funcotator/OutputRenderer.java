package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.LinkedHashMap;
import java.util.List;

/**
 * An abstract class to allow for writing output for the Funcotator.
 * Used to output Funcotations to a location of the user's choice.
 * For example, writing out to a VCF file (e.g. {@link org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer}).
 * It is also possible to write to a stream, or other connection (e.g. a UDP port)
 * at the user's discretion.
 * Created by jonn on 8/30/17.
 */
public abstract class OutputRenderer implements AutoCloseable {

    //==================================================================================================================
    /**
     * {@link LinkedHashMap} of manually specified annotations to add to each output in addition to annotations provided
     * to {@link OutputRenderer#write(VariantContext, List)}.
     */
    protected LinkedHashMap<String, String> manualAnnotations;

    /**
     * {@link String} representation of {@link OutputRenderer#manualAnnotations} serialized to the output format of this {@link OutputRenderer}.
     */
    protected String manualAnnotationSerializedString;

    //==================================================================================================================

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
