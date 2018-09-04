package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.ojalgo.commons.math3.linear.Access2DWrapper;
import org.ojalgo.commons.math3.linear.RealMatrixWrapper;
import org.ojalgo.matrix.decomposition.SingularValue;

/**
 * SVD using the ojAlgo library.
 */
public final class OjAlgoSingularValueDecomposer implements SingularValueDecomposer {

    private static final Logger logger = LogManager.getLogger(OjAlgoSingularValueDecomposer.class);

    /**
     * Create a SVD instance using ojAlgo.
     *
     * @param m matrix that is not {@code null}
     * @return SVD instance that is never {@code null}
     */
    @Override
    public SVD createSVD(final RealMatrix m) {

        Utils.nonNull(m, "Cannot create SVD on a null matrix.");

        logger.info("Converting Apache RealMatrix to ojAlgo Matrix...");
        final RealMatrixWrapper pds = RealMatrixWrapper.of(m);

        logger.info("Calculating SVD...");
        final SingularValue<Double> svd = SingularValue.make(pds);
        svd.compute(pds);

        logger.info("Converting ojAlgo Matrix to Apache RealMatrix...");
        final RealMatrix pinv = Access2DWrapper.of(svd.getInverse());
        final RealMatrix u = Access2DWrapper.of(svd.getQ1());
        final RealMatrix v = Access2DWrapper.of(svd.getQ2());
        final double[] s = svd.getSingularValues().toRawCopy1D();

        return new SimpleSVD(u, s, v, pinv);
    }
}
