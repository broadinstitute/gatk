package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.List;

/**
 * An abstract class to allow for the creation of a {@link Funcotation} for a given data source.
 * Created by jonn on 8/30/17.
 */
public interface DataSourceFuncotationFactory extends AutoCloseable {

    /**
     * Perform cleanup tasks for this {@link DataSourceFuncotationFactory}.
     */
    default void close() {}

    /**
     * @return An ordered list of the names of annotations that this Data Source supports.
     */
    List<String> getSupportedFuncotationFields();

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, and {@code featureContext}.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureList {@link List} of {@link Feature} corresponding to the given {@code variant}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList);
}
