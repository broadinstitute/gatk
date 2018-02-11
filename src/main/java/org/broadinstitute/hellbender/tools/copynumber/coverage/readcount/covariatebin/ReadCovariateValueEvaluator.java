package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
@FunctionalInterface
public interface ReadCovariateValueEvaluator {

    /**
     * Get value of the covariate for a given read, and reference context surrounding it
     *
     * @param read GATK read
     * @param referenceContext reference context surrounding the read
     * @return covariate value
     */
    double getCovariateValue(final GATKRead read, final ReferenceContext referenceContext);
}
