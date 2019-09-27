package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;

import java.util.*;

public class ReferenceContextFuncotationFactory extends ComputedFuncotationFactory {

    private static final long serialVersionUID = 11L;

    //==================================================================================================================
    // Public Static Members:
    /**
     * Default name for this funcotation data source.
     */
    public static final String DEFAULT_NAME = "ComputedReferenceContext";

    /**
     * Default version for this funcotation data source.
     */
    public static final String DEFAULT_VERSION = "1";

    /**
     * The field name for the reference context annotation
     */
    public static final String REFERENCE_CONTEXT_FIELD_NAME = "ReferenceContext";

    //==================================================================================================================
    // Private Static Members:
    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(ReferenceContextFuncotationFactory.class);

    /**
     * Name for this funcotation data source.
     */
    private final String name;

    /**
     * The window around a variant to include in the reference context annotation.
     */
    private final int referenceContextWindowSize;


    /**
     * Creates a {@link ReferenceContextFuncotationFactory}.
     *
     * @param version             Version of this {@link ReferenceContextFuncotationFactory}.
     * @param name                A {@link String} containing the name of this {@link ReferenceContextFuncotationFactory}.
     * @param referenceContextWindowSize     The window around a variant to include in the reference context annotation.
     */
    public ReferenceContextFuncotationFactory(String version, String name, int referenceContextWindowSize) {
        this.version = version;
        this.name = name;
        this.referenceContextWindowSize = referenceContextWindowSize;
    }


    @Override
    public String getName() {
        return name;
    }

    /**
     * Creates a {@link List} of {@link Funcotation}s for the given {@code variant}, {@code referenceContext}, {@code featureContext}, and {@code gencodeFuncotations}.
     * For some Data Sources knowledge of Gene Name or Transcript ID is required for annotation.
     * Accounts for override values passed into the constructor as well.
     * @param variant {@link VariantContext} to annotate.  Never {@code null}.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.  Never {@code null}.
     * @param featureContext {@link FeatureContext} corresponding to the variant.  Never {@code null}.
     * @param gencodeFuncotations {@link List} of {@link GencodeFuncotation} that have already been created for the given {@code variant}/{@code referenceContext}/{@code featureContext}.
     *   {@code null} is acceptable since there is no dependency on Gencode annotations.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    @Override
    public List<Funcotation> createFuncotations(VariantContext variant, ReferenceContext referenceContext, FeatureContext featureContext, List<GencodeFuncotation> gencodeFuncotations) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        for ( final Allele altAllele : variant.getAlternateAlleles()) {

            // Get the reference bases for our current variant:
            final StrandCorrectedReferenceBases referenceBases = FuncotatorUtils.createReferenceSnippet(variant.getReference(), altAllele, referenceContext, Strand.POSITIVE, referenceContextWindowSize);

            TableFuncotation referenceContextFuncotation = TableFuncotation.create(
                    Collections.singletonList(REFERENCE_CONTEXT_FIELD_NAME),
                    Collections.singletonList(referenceBases.getBaseString(Strand.POSITIVE)),
                    altAllele,
                    getName(),
                    FuncotationMetadataUtils.createWithUnknownAttributes(Collections.singletonList(REFERENCE_CONTEXT_FIELD_NAME)));
            outputFuncotations.add(referenceContextFuncotation);
        }

        return outputFuncotations;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        LinkedHashSet<String> supportedFuncotationFields = new LinkedHashSet<>();
        supportedFuncotationFields.add(REFERENCE_CONTEXT_FIELD_NAME);
        return supportedFuncotationFields;
    }

    @Override
    public Class<? extends Feature> getAnnotationFeatureClass() {
        // Returning Feature.class here implies that this class doesn't care about what features it gets.
        return Feature.class;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFieldsForSegments() {
        return new LinkedHashSet<String>();
    }

    @Override
    public boolean isSupportingSegmentFuncotation() {
        return false;
    }
}
