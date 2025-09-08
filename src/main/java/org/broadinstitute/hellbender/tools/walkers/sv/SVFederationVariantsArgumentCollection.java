package org.broadinstitute.hellbender.tools.walkers.sv;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.util.Arrays;
import java.util.List;

public class SVFederationVariantsArgumentCollection extends MultiVariantInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    public static final String VARIANTS_A_LONG_NAME = "variants-A";
    public static final String VARIANTS_A_SHORT_NAME = "A";
    @Argument(
            doc = "Variants file from cohort A in VCF format",
            fullName = VARIANTS_A_LONG_NAME,
            shortName = VARIANTS_A_SHORT_NAME
    )
    public FeatureInput<VariantContext> varA;

    public static final String VARIANTS_B_LONG_NAME = "variants-B";
    public static final String VARIANTS_B_SHORT_NAME = "B";
    @Argument(
            doc = "Variants file from cohort B in VCF format",
            fullName = VARIANTS_B_LONG_NAME,
            shortName = VARIANTS_B_SHORT_NAME
    )
    public FeatureInput<VariantContext> varB;


    @Override
    public List<GATKPath> getDrivingVariantPaths() {
        return Arrays.asList(varA, varB);
    }

    public List<FeatureInput<VariantContext>> getFeatureInputsForDrivingVariants() { return Arrays.asList(varA, varB); }
}
