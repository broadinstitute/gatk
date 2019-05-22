package org.broadinstitute.hellbender.tools.funcotator.compositeoutput;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.funcotator.FuncotationMap;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;

import java.util.List;

/**
 *  Class to make multiple funcotator output at the same time.
 */
public class CompositeOutputRenderer extends OutputRenderer {

    private final List<OutputRenderer> outputRenderers;

    public CompositeOutputRenderer(final List<OutputRenderer> outputRenderers, final String toolVersion) {
        super(toolVersion);
        this.outputRenderers = outputRenderers;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() {
        outputRenderers.forEach(OutputRenderer::close);
    }

    /**
     * {@inheritDoc}
     *
     * @param variant {@link VariantContext} to write to the file.
     * @param txToFuncotationMap {@link FuncotationMap} to add to the given {@code variant} on output.
     */
    @Override
    public void write(final VariantContext variant, final FuncotationMap txToFuncotationMap) {
        outputRenderers.forEach(or -> or.write(variant, txToFuncotationMap));
    }
}
