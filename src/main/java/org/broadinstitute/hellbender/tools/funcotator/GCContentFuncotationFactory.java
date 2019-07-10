package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadataUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;

public class GCContentFuncotationFactory extends ComputedFuncotationFactory {

    private static final long serialVersionUID = 11L;

    //==================================================================================================================
    // Public Static Members:
    /**
     * Default name for this funcotation data source.
     */
    public static final String DEFAULT_NAME = "ComputedGCContent";

    /**
     * Default version for this funcotation data source.
     */
    public static final String DEFAULT_VERSION = "1";

    /**
     * The field name for the GC content annotation
     */
    public static final String GC_CONTENT_FIELD_NAME = "GCContent";

    /**
     * Name for this funcotation data source.
     */
    private final String name;

    /**
     * Number of bases to the left and right of a variant in which to calculate the GC content.
     */
    private final int gcContentWindowSize;


    /**
     * Creates a {@link GCContentFuncotationFactory}.
     *
     * @param version             Version of this {@link GCContentFuncotationFactory}.
     * @param name                A {@link String} containing the name of this {@link GCContentFuncotationFactory}.
     * @param gcContentWindowSize Number of bases around the locus of the variant to calculate the GC content.
     */
    public GCContentFuncotationFactory(String version, String name, int gcContentWindowSize) {
        this.version = version;
        this.name = name;
        this.gcContentWindowSize = gcContentWindowSize;
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
     * @param gencodeFuncotations {@link List} of {@link GencodeFuncotation}s that have already been created for the given {@code variant}/{@code referenceContext}/{@code featureContext}.
     *   {@code null} is acceptable since there is no dependency on Gencode annotations.
     * @return {@link List} of {@link Funcotation}s given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    @Override
    public List<Funcotation> createFuncotations(VariantContext variant, ReferenceContext referenceContext, FeatureContext featureContext, List<GencodeFuncotation> gencodeFuncotations) {
        final List<Funcotation> outputFuncotations = new ArrayList<>();

        for ( final Allele altAllele : variant.getAlternateAlleles()) {
            TableFuncotation gcContentFuncotation = TableFuncotation.create(
                    Collections.singletonList(GC_CONTENT_FIELD_NAME),
                    Collections.singletonList(String.valueOf(calculateGcContent(variant.getReference(), altAllele, referenceContext, gcContentWindowSize))),
                    altAllele,
                    getName(),
                    FuncotationMetadataUtils.createWithUnknownAttributes(Collections.singletonList(GC_CONTENT_FIELD_NAME)));
            outputFuncotations.add(gcContentFuncotation);
        }

        return outputFuncotations;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        LinkedHashSet<String> supportedFuncotationFields = new LinkedHashSet<>();
        supportedFuncotationFields.add(GC_CONTENT_FIELD_NAME);
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


    /**
     * Calculates the fraction of Guanine and Cytosine bases in a window of a given size around a variant.
     * Note: Since Guanine and Cytosine are complementary bases, strandedness makes no difference.
     * @param refAllele Reference {@link Allele} for the locus in question.
     * @param altAllele Alternate {@link Allele} for the locus in question.
     * @param referenceContext The {@link ReferenceContext} for a variant.  Assumed to already be centered on the variant of interest.  Must not be {@code null}.
     * @param windowSize The number of bases to the left and right of the given {@code variant} to calculate the GC Content.  Must be >=1.
     * @return The fraction of Guanine and Cytosine bases / total bases in a window of size {@code windowSize} around a variant.
     */
    public static double calculateGcContent( final Allele refAllele,
                                             final Allele altAllele,
                                             final ReferenceContext referenceContext,
                                             final int windowSize ) {

        // TODO: this seems to do something similar to FuncotatorUtils::getBasesInWindowAroundReferenceAllele - should this method call into that?

        Utils.nonNull( referenceContext );
        ParamUtils.isPositive( windowSize, "Window size must be >= 1." );

        final int leadingWindowSize;
        final int trailingWindowSize = windowSize;

        if ( GATKVariantContextUtils.isInsertion(refAllele, altAllele) ||
                GATKVariantContextUtils.isDeletion(refAllele, altAllele)) {
            // If we have an insertion, we take 1 less base from the front
            // because the insertion happens between two codons.
            // The preceding padding base is there as a convenience in VCF files.
            // Thus the prior <windowSize> bases will contain this leading padding base.

            // If we have a deletion, the convention in VCF files is to include a
            // padding base at the front prior to the deleted bases so the alternate
            // allele can be non-empty.
            // Because of this we subtract 1 from the leading window size.

            leadingWindowSize = windowSize - 1;
        }
        else {
            leadingWindowSize = windowSize;
        }

        // Get the bases:
        byte[] bases = referenceContext.getBases(leadingWindowSize, trailingWindowSize);

        // Get the gcCount:
        long gcCount = 0;
        for ( final byte base : bases ) {
            if ( BaseUtils.basesAreEqual(base, BaseUtils.Base.G.base) || BaseUtils.basesAreEqual(base, BaseUtils.Base.C.base) ) {
                ++gcCount;
            }
        }

        // Calculate the ratio itself:
        return ((double)gcCount) / ((double) bases.length);
    }
}
