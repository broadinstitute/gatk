package org.broadinstitute.hellbender.tools.walkers.sv;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import javax.ws.rs.core.Variant;
import java.util.Arrays;
import java.util.List;

public class SVFederationArgumentCollection extends MultiVariantInputArgumentCollection {
    private static final long serialVersionUID = 1L;

    /**
     * Expected format is tab-delimited and contains columns VID_A, VID_B, SCORE.
     * First line must be a header with column names. Comment lines starting with
     * {@link TableUtils#COMMENT_PREFIX} are ignored.
     */
    public static final String SV_PAIR_FILE_LONG_NAME = "sv-pairs";
    @Argument(
            doc = "SV pair file (.tsv) containing the candidate SV pairs and matching scores",
            fullName = SV_PAIR_FILE_LONG_NAME
    )
    public GATKPath svPairFilePath;

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

    public final GATKPath getSVPairFilePath() { return svPairFilePath; }
}
