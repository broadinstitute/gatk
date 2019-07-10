package org.broadinstitute.hellbender.tools.funcotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.util.LinkedHashSet;
import java.util.List;

public abstract class FuncotationFactory {

    /** Default version string for this {@link FuncotationFactory}. */
    @VisibleForTesting
    public static final String DEFAULT_VERSION_STRING = "UNKNOWN_VERSION";

    /**
     * Version number of this {@link FuncotationFactory}.
     */
    protected String version = DEFAULT_VERSION_STRING;

    /**
     * @return A {@link String} containing information about this {@link FuncotationFactory}.
     */
    public String getInfoString() {
        return getName() + " " + getVersion();
    }

    /**
     * @return The name of the data source corresponding to this {@link FuncotationFactory}.
     */
    public abstract String getName();

    /**
     * @return The version of the data source corresponding to this {@link FuncotationFactory}.
     */
    public String getVersion() {
        return version;
    }

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
    public abstract List<Funcotation> createFuncotations(final VariantContext variant, final ReferenceContext referenceContext, final FeatureContext featureContext, final List<GencodeFuncotation> gencodeFuncotations);


    /**
     * @return Get the {@link Class} of the feature type that can be used to create annotations by this {@link FuncotationFactory}.
     */
    @VisibleForTesting
    public abstract Class<? extends Feature> getAnnotationFeatureClass();

    /**
     * @return An ordered {@link LinkedHashSet} of the names of annotations that this Data Source supports.
     */
    public abstract LinkedHashSet<String> getSupportedFuncotationFields();

    /**
     * @return An ordered {@link LinkedHashSet} of the names of annotations that this Data Source supports when annotating segments.
     */
    public abstract LinkedHashSet<String> getSupportedFuncotationFieldsForSegments();


    /**
     * @return Whether this funcotation factory can support creating funcotations from segments.
     */
    public abstract boolean isSupportingSegmentFuncotation();


}
