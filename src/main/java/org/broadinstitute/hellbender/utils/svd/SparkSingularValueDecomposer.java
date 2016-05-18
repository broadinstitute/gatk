package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.*;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.spark.SparkConverter;

import java.util.Arrays;

class SparkSingularValueDecomposer {

    private static final double EPS = 1e-32;

    private static final Logger logger = LogManager.getLogger(SparkSingularValueDecomposer.class);

    private final static int NUM_SLICES = 60;

    /**
     *  Create a SVD of the given matrix using the given Java Spark Context.
     *
     * @param sc SparkContext to run the SVD. Not {@code null}
     * @param realMat the matrix target.  Not {@code null}
     * @return never {@code null}
     */
    public static SVD createSVD(final JavaSparkContext sc, final RealMatrix realMat){

        Utils.nonNull(sc, "Cannot perform Spark MLLib SVD using a null JavaSparkContext.");
        Utils.nonNull(realMat, "Cannot perform Spark MLLib SVD on a null matrix.");

        final RowMatrix mat = SparkConverter.convertRealMatrixToSparkRowMatrix(sc, realMat, NUM_SLICES);

        // Compute all of the singular values and corresponding singular vectors.
        final SingularValueDecomposition<RowMatrix, Matrix> svd = mat.computeSVD((int) mat.numCols(), true, 1.0E-9d);

        // Get our distributed results
        final RowMatrix u = svd.U();
        final Vector s = svd.s();
        final Matrix v = svd.V().transpose();

        // Move the matrices from Spark/distributed space to Apache Commons space
        logger.info("Converting distributed Spark matrix to local matrix...");
        final RealMatrix uReal = SparkConverter.convertSparkRowMatrixToRealMatrix(u, realMat.getRowDimension());
        logger.info("Done converting distributed Spark matrix to local matrix...");
        logger.info("Converting Spark matrix to local matrix...");
        final RealMatrix vReal = SparkConverter.convertSparkMatrixToRealMatrix(v);
        logger.info("Done converting Spark matrix to local matrix...");
        final double [] singularValues = s.toArray();

        logger.info("Calculating the pseudoinverse...");
        logger.info("Pinv: calculating tolerance...");

        // Note that the pinv of realMat is V * invS * U'
        final double tolerance = Math.max(realMat.getColumnDimension(), realMat.getRowDimension()) * realMat.getNorm() * EPS;
        logger.info("Pinv: inverting the singular values (with tolerance) and creating a diagonal matrix...");
        final double[] invS = Arrays.stream(singularValues).map(sv -> invertSVWithTolerance(sv, tolerance)).toArray();

        final Matrix invSMat = Matrices.diag(Vectors.dense(invS));
        logger.info("Pinv: Multiplying V * invS * U' to get the pinv (using pinv transpose = U * invS' * V') ...");
        final RowMatrix pinvT = u.multiply(invSMat).multiply(v);
        logger.info("Pinv: Converting back to local matrix ...");
        final RealMatrix pinv = SparkConverter.convertSparkRowMatrixToRealMatrix(pinvT, realMat.getRowDimension()).transpose();
        logger.info("Done calculating the pseudoinverse and converting it...");

        return new SimpleSVD(uReal, s.toArray(), vReal, pinv);
    }

    private static double invertSVWithTolerance(final double sv, final double tol){
        if (sv <= tol) {
            return 0;
        } else {
            return 1/sv;
        }
    }
}
