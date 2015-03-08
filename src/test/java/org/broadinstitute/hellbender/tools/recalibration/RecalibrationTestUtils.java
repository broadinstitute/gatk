package org.broadinstitute.hellbender.tools.recalibration;

import org.broadinstitute.hellbender.tools.recalibration.covariates.*;

import java.util.Arrays;
import java.util.List;

public class RecalibrationTestUtils {
    public static List<Covariate> makeInitializedStandardCovariates() {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
        final List<Covariate> covariates = Arrays.asList(new ReadGroupCovariate(), new QualityScoreCovariate(), new ContextCovariate(), new CycleCovariate());
        for ( Covariate cov : covariates ) {
            cov.initialize(RAC);
        }
        return covariates;
    }
}
