package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.io.Closeable;
import java.util.*;

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
     * The backing data store as a FeatureInput to leverage tribble querying.  Can be {@code null} for non-locatable
     * funcotation factories.
     */
    protected final FeatureInput<? extends Feature> mainSourceFileAsFeatureInput;

    @VisibleForTesting
    public FeatureInput<? extends Feature> getMainSourceFileAsFeatureInput() {
        return mainSourceFileAsFeatureInput;
    }

    /**
     * Constructor to initialize final fields in this class with defaults.
     */
    protected DataSourceFuncotationFactory() {
        this.mainSourceFileAsFeatureInput = null;
    }

    /**
     * Constructor to initialize final fields in this class.
     * @param mainSourceFileAsFeatureInput The backing data store as a FeatureInput to leverage tribble querying.  Can be {@code null} for non-locatable funcotation factories.
     */
    protected DataSourceFuncotationFactory(final FeatureInput<? extends Feature> mainSourceFileAsFeatureInput) {
        this.mainSourceFileAsFeatureInput = mainSourceFileAsFeatureInput;
    }


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
     * @param featureContext {@link FeatureContext} corresponding to the variant.  Never {@code null}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        return createFuncotations(variant, referenceContext, featureContext, null);
    }

    /**
     * Creates a {@link List} of {@link Funcotation} for the given {@code variant}, {@code referenceContext}, {@code featureContext}, and {@code gencodeFuncotations}.
     * For some Data Sources knowledge of Gene Name or Transcript ID is required for annotation.
     * Accounts for override values passed into the constructor as well.
     * @param variant {@link VariantContext} to annotate.  Never {@code null}.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.  Never {@code null}.
     * @param featureContext {@link FeatureContext} corresponding to the variant.  Never {@code null}.
     * @param gencodeFuncotations {@link List} of {@link GencodeFuncotation} that have already been created for the given {@code variant}/{@code referenceContext}/{@code featureContext}.
     *   {@code null} is acceptable if there are no corresponding gencode funcotations.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    public List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext, final List<GencodeFuncotation> gencodeFuncotations) {

        Utils.nonNull(variant);
        Utils.nonNull(referenceContext);
        Utils.nonNull(featureContext);

        // Query this funcotation factory to get the list of overlapping features.
        final List<Feature> featureList = queryFeaturesFromFeatureContext(featureContext);

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
     * Queries the provided FeatureContext for Features from our FeatureInput {@link #mainSourceFileAsFeatureInput}.
     * The default implementation returns all Features from our FeatureInput that overlap the FeatureContext's
     * interval, but subclasses may override (for example, to pad the query).
     *
     * @param featureContext the FeatureContext to query
     * @return Features from our FeatureInput {@link #mainSourceFileAsFeatureInput} queried from the FeatureContext
     */
    protected List<Feature> queryFeaturesFromFeatureContext(final FeatureContext featureContext) {
        @SuppressWarnings("unchecked")
        final List<Feature> features = (List<Feature>)featureContext.getValues(mainSourceFileAsFeatureInput);
        return features;
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
    @VisibleForTesting
    public abstract Class<? extends Feature> getAnnotationFeatureClass();
}
