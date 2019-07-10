package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotationFactory;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.gencode.GencodeGtfGeneFeature;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.file.Path;
import java.util.*;

public abstract class ComputedFuncotationFactory extends FuncotationFactory {

    //==================================================================================================================
    // Public Static Members:
    /**
     * Default name for this data source (i.e. computed in this case).
     */
    public static final String DEFAULT_NAME = "ComputedFuncotations";


    /**
     * Creates a {@link List} of {@link Funcotation}s for the given {@code variant}, {@code referenceContext}, {@code featureContext}, and {@code gencodeFuncotations}.
     * For some Data Sources knowledge of Gene Name or Transcript ID is required for annotation.
     * Accounts for override values passed into the constructor as well.
     * @param variant {@link VariantContext} to annotate.  Never {@code null}.
     * @param referenceContext {@link ReferenceContext} corresponding to the given {@code variant}.  Never {@code null}.
     * @param featureContext {@link FeatureContext} corresponding to the variant.  Never {@code null}.
     * @return {@link List} of {@link Funcotation} given the {@code variant}, {@code referenceContext}, and {@code featureContext}.  This should never be empty.
     */
    @Override
    public List<Funcotation> createFuncotations(VariantContext variant, ReferenceContext referenceContext, FeatureContext featureContext) {
        return createFuncotations(variant, referenceContext, featureContext, null);
    }
}
