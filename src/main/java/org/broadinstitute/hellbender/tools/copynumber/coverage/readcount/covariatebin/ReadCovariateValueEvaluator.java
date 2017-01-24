package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.covariatebin;

import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
@FunctionalInterface
public interface ReadCovariateValueEvaluator {

    /**
     * Get value of the covariate for a given read
     *
     * @param read GATK read
     * @return covariate value
     */
    double getValueFromRead(final GATKRead read);
}
