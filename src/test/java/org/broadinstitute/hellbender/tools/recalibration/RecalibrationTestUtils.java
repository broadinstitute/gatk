package org.broadinstitute.hellbender.tools.recalibration;

import org.broadinstitute.hellbender.tools.recalibration.covariates.*;

public final class RecalibrationTestUtils {
    public static Covariate[] makeInitializedStandardCovariates() {
        final RecalibrationArgumentCollection RAC = new RecalibrationArgumentCollection();
        final Covariate[] covariates = new Covariate[4];
        covariates[0] = new ReadGroupCovariate();
        covariates[1] = new QualityScoreCovariate();
        covariates[2] = new ContextCovariate();
        covariates[3] = new CycleCovariate();
        for ( Covariate cov : covariates ) cov.initialize(RAC);
        return covariates;
    }
}
