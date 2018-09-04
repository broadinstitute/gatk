package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.io.Closeable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * An abstract class to allow for the creation of a {@link Funcotation} for a given data source.
 * Created by jonn on 8/30/17.
 */
public abstract class DataSourceFuncotationFactory implements Closeable {

    //==================================================================================================================

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(DataSourceFuncotationFactory.class);

    /** Default version string for this {@link DataSourceFuncotationFactory}. */
    @VisibleForTesting
    public static final String DEFAULT_VERSION_STRING = "UNKNOWN_VERSION";

    /**
     * Version number of this {@link DataSourceFuncotationFactory}.
     */
    protected String version = DEFAULT_VERSION_STRING;

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
     * @return A {@link String} containing information about this {@link DataSourceFuncotationFactory}.
     */
    public String getInfoString() {
        return getName() + " " + getVersion();
    }

    /**
     * Perform cleanup tasks for this {@link DataSourceFuncotationFactory}.
     */
    public void close() {}

    /**
     * Apply the override values in {@link DataSourceFuncotationFactory#annotationOverrideMap} to every
     * {@link Funcotation} in the given {@code outputFuncotations}.
     * @param funcotationList {@link List} of {@link Funcotation} to which to apply override values.
     */
    protected void setOverrideValuesInFuncotations(final List<Funcotation> funcotationList) {
        for ( final Funcotation funcotation : funcotationList ) {
            funcotation.setFieldSerializationOverrideValues( annotationOverrideMap );
        }
    }

    /**
     * @return The name of the data source corresponding to this {@link DataSourceFuncotationFactory}.
     */
    public abstract String getName();

    /**
     * @return The {@link org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions.DataSourceType} of this {@link DataSourceFuncotationFactory}.
     */
    public abstract FuncotatorArgumentDefinitions.DataSourceType getType();

    /**
     * @return The version of the data source corresponding to this {@link DataSourceFuncotationFactory}.
     */
    public String getVersion() {
        return version;
    }

    /**
     * @return An ordered {@link LinkedHashSet} of the names of annotations that this Data Source supports.
     */
    public abstract LinkedHashSet<String> getSupportedFuncotationFields();

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, and {@code featureContext}.
     * Accounts for override values passed into the constructor as well.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureSourceMap {@link Map} of {@link String} -> {@link List} of {@link Feature} (data source name -> data source features corresponding to the given {@code variant}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final Map<String, List<Feature>> featureSourceMap) {
        return createFuncotations(variant, referenceContext, featureSourceMap, null);
    }

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, {@code featureContext}, and {@code gencodeFuncotations}.
     * For some Data Sources knowledge of Gene Name or Transcript ID is required for annotation.
     * Accounts for override values passed into the constructor as well.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureSourceMap {@link Map} of {@link String} -> {@link List} of {@link Feature} (data source name -> data source features corresponding to the given {@code variant}.
     * @param gencodeFuncotations {@link List} of {@link GencodeFuncotation} that have already been created for the given {@code variant}/{@code referenceContext}/{@code featureContext}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final Map<String, List<Feature>> featureSourceMap, final List<GencodeFuncotation> gencodeFuncotations) {

        // Get the features that this funcotation factory is responsible for:
        final List<Feature> featureList = getFeatureListFromMap(featureSourceMap);

        final List<Funcotation> outputFuncotations;

        // If our featureList is compatible with this DataSourceFuncotationFactory, then we make our funcotations:
        if ( isFeatureListCompatible(featureList) ) {

            // Create our funcotations:
            if ( gencodeFuncotations == null ) {
                outputFuncotations = createFuncotationsOnVariant(variant, referenceContext, featureList);
            }
            else {
                outputFuncotations = createFuncotationsOnVariant(variant, referenceContext, featureList, gencodeFuncotations);
            }

            // Set our overrides:
            setOverrideValuesInFuncotations(outputFuncotations);
        }
        else {
            return createDefaultFuncotationsOnVariant(variant, referenceContext);
        }

        if ((outputFuncotations == null) || (outputFuncotations.size() == 0)) {
            return createDefaultFuncotationsOnVariant(variant, referenceContext);
        } else {
            return outputFuncotations;
        }
    }

    /**
     * Checks to see if the given featureList is compatible with this {@link DataSourceFuncotationFactory}.
     * Cues off of the feature type in the feature list and whether the given list contains any non-null features.
     * This method acts as a sanity-check before attempting to do any annotations on features.
     * @param featureList {@link List} of {@link Feature} that might be applicable to this {@link DataSourceFuncotationFactory} for annotation.
     * @return {@code true} if the given {@code featureList} contains at least one non-null feature of type {@link #getAnnotationFeatureClass()}; {@code false} otherwise.
     */
    private boolean isFeatureListCompatible(final List<Feature> featureList) {
        // Make sure these features can be annotated by this DataSourceFuncotationFactory:
        // NOTE: We only check the first non-null element of the list for feature type:
        boolean foundCompatibleFeature = false;
        for ( final Feature f : featureList ) {
            if (f != null) {
                foundCompatibleFeature = getAnnotationFeatureClass().isAssignableFrom(f.getClass());
                break;
            }
        }
        return foundCompatibleFeature;
    }

    /**
     * Get the list of features to annotate from the given Map of features.
     * Extracts the feature list given the name of this {@link DataSourceFuncotationFactory}.
     * @param featureSourceMap {@link Map} of {@link String} -> ({@link List} of {@link Feature}) (Data source name -> feature list) containing all features that could be used for this {@link DataSourceFuncotationFactory}.
     * @return A {@link List} of {@link Feature} that are to be annotated by this {@link DataSourceFuncotationFactory}
     */
    private List<Feature> getFeatureListFromMap(final Map<String, List<Feature>> featureSourceMap) {
        // Get the features that this funcotation factory is responsible for:
        final List<Feature> featureList;

        // Only worry about name filtering if we care about the specific feature type:
        // NOTE: This should probably be fixed to key off some other abstract class logic.
        if ( getAnnotationFeatureClass().equals(Feature.class) ) {
            featureList = featureSourceMap.entrySet().stream()
                    .map(Map.Entry::getValue)
                    .flatMap(Collection::stream)
                    .collect(Collectors.toList());
        }
        else {
            featureList = featureSourceMap.getOrDefault( getName(), Collections.emptyList() );
        }
        return featureList;
    }

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant} and {@code referenceContext}.
     * These will be default funcotations that essentially have empty values.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    protected abstract List<Funcotation> createDefaultFuncotationsOnVariant( final VariantContext variant, final ReferenceContext referenceContext);

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, and {@code featureContext}.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureList {@link List} of {@link Feature} corresponding to the given {@code variant}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    protected abstract List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList);

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, {@code featureContext}, and {@code gencodeFuncotations}.
     * For some Data Sources knowledge of Gene Name or Transcript ID is required for annotation.
     * @param variant {@link VariantContext} to annotate.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.
     * @param featureList {@link List} of {@link Feature} corresponding to the given {@code variant}.
     * @param gencodeFuncotations {@link List} of {@link GencodeFuncotation} that have already been created for the given {@code variant}/{@code referenceContext}/{@code featureContext}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    protected abstract List<Funcotation> createFuncotationsOnVariant(final VariantContext variant,
                                                                  final ReferenceContext referenceContext,
                                                                  final List<Feature> featureList,
                                                                  final List<GencodeFuncotation> gencodeFuncotations);

    /**
     * @return Get the {@link Class} of the feature type that can be used to create annotations by this {@link DataSourceFuncotationFactory}.
     */
    protected abstract Class<? extends Feature> getAnnotationFeatureClass();
}
