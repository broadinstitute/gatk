package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;

import java.util.*;

/**
 * An abstract class to allow for the creation of a {@link Funcotation} for a given data source.
 * Created by jonn on 8/30/17.
 */
public abstract class DataSourceFuncotationFactory implements AutoCloseable {

    //==================================================================================================================

    /**
     * Map of ANNOTATION_NAME -> OVERRIDE_VALUE.
     */
    protected Map<String, String> annotationOverrideMap;

    /**
     * Set values in {@link DataSourceFuncotationFactory#annotationOverrideMap} based on the given annotation override values
     * and whether or not this {@link DataSourceFuncotationFactory} supports those annotations.
     * @param annotationOverrides The {@link Map} of annotation override key names and values.
     */
    protected void initializeAnnotationOverrides(final LinkedHashMap<String, String> annotationOverrides) {
        // Go through the Annotation Maps and check to see if the default/override annotation names are applicable for
        // this FuncotationFactory:
        final LinkedHashSet<String> supportedFuncotations = new LinkedHashSet<>( getSupportedFuncotationFields() );
        this.annotationOverrideMap = new HashMap<>();
        for ( final String annotationOverrideKey : annotationOverrides.keySet() ) {
            if ( supportedFuncotations.contains(annotationOverrideKey) ) {
                annotationOverrideMap.put( annotationOverrideKey, annotationOverrides.get(annotationOverrideKey) );
            }
        }
    }

    //==================================================================================================================

    /**
     * Perform cleanup tasks for this {@link DataSourceFuncotationFactory}.
     */
    public void close() {}

    /**
     * @return The name of the data source corresponding to this {@link DataSourceFuncotationFactory}.
     */
    public abstract String getName();

    /**
     * @return An ordered {@link LinkedHashSet} of the names of annotations that this Data Source supports.
     */
    public abstract LinkedHashSet<String> getSupportedFuncotationFields();

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, and {@code featureContext}.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureList {@link List} of {@link Feature} corresponding to the given {@code variant}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    public abstract List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList);
}
