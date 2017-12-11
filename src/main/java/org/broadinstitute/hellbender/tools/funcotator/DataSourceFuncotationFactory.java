package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

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
        final Set<String> supportedFuncotations = getSupportedFuncotationFields();
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
     * Apply the override values in {@link DataSourceFuncotationFactory#annotationOverrideMap} to every
     * {@link Funcotation} in the given {@code outputFuncotations}.
     * @param funcotationList {@link List} of {@link Funcotation} to which to apply override values.
     */
    public void setOverrideValuesInFuncotations(final List<Funcotation> funcotationList) {
        for ( final Funcotation funcotation : funcotationList ) {
            funcotation.setFieldSerializationOverrideValues( annotationOverrideMap );
        }
    }

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

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, {@code featureContext}, and {@code gencodeFuncotations}.
     * For some Data Sources knowledge of Gene Name or Transcript ID is required for annotation.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureList {@link List} of {@link Feature} corresponding to the given {@code variant}.
     * @param gencodeFuncotations {@link List} of {@link GencodeFuncotation} that have already been created for the given {@code variant}/{@code referenceContext}/{@code featureContext}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    public abstract List<Funcotation> createFuncotations(final VariantContext variant,
                                                         final ReferenceContext referenceContext,
                                                         final List<Feature> featureList,
                                                         final List<GencodeFuncotation> gencodeFuncotations);
}
