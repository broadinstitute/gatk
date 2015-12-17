package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.*;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.LinkedList;


class SparkSingularValueDecomposer {

    private static final double EPS = 1e-32;

    private static final Logger logger = LogManager.getLogger(SparkSingularValueDecomposer.class);

    private final static int NUM_SLICES = 50;

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

        logger.info("Converting matrix to distributed Spark matrix...");
        final double [][] dataArray = realMat.getData();
        final LinkedList<Vector> rowsList = new LinkedList<>();
        for (final double [] i : dataArray) {
            final Vector currentRow = Vectors.dense(i);
            rowsList.add(currentRow);
        }

        // We may want to swap out this static value for something dynamic (as shown below), but this seems to slow it down.
        // final int totalSpace = realMat.getColumnDimension() * realMat.getRowDimension() * Double.BYTES;
        // // Want the partitions to be ~100KB of space
        // final int slices = totalSpace/100000;
        final JavaRDD<Vector> rows = sc.parallelize(rowsList, NUM_SLICES);

        // Create a RowMatrix from JavaRDD<Vector>.
        final RowMatrix mat = new RowMatrix(rows.rdd());
        logger.info("Done converting matrix to distributed Spark matrix...");

        // Compute all of the singular values and corresponding singular vectors.
        final SingularValueDecomposition<RowMatrix, Matrix> svd = mat.computeSVD((int) mat.numCols(), true, 1.0E-9d);

        // Get our distributed results
        final RowMatrix u = svd.U();
        final Vector s = svd.s();
        final Matrix v = svd.V().transpose();

        // Move the matrices from Spark/distributed space to Apache Commons space
        logger.info("Converting distributed Spark matrix to local matrix...");
        final RealMatrix uReal = convertSparkRowMatrixToRealMatrix(u, realMat.getRowDimension());
        logger.info("Done converting distributed Spark matrix to local matrix...");
        logger.info("Converting Spark matrix to local matrix...");
        final RealMatrix vReal = convertSparkMatrixToRealMatrix(v);
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
        final RowMatrix pinvT = (u.multiply(invSMat)).multiply(v);
        final RealMatrix pinv = convertSparkRowMatrixToRealMatrix(pinvT, realMat.getRowDimension()).transpose();
        logger.info("Done calculating the pseudoinverse and converting it...");

        return new SimpleSVD(uReal, s.toArray(), vReal, pinv);
    }

    private static RealMatrix convertSparkMatrixToRealMatrix(final Matrix r){
        final RealMatrix result = new Array2DRowRealMatrix(r.numRows(), r.numCols());
        final double [] columnMajorMat = r.toArray();
        for (int i = 0; i < r.numRows(); i++) {
            result.setRow(i, Arrays.copyOfRange(columnMajorMat, i * r.numCols(), i * r.numCols() + r.numCols()) );
        }
        return result;
    }

    private static RealMatrix convertSparkRowMatrixToRealMatrix(final RowMatrix r, final int cachedNumRows) {

        if (r == null) {
            return null;
        }

        int numRows;
        if (cachedNumRows == -1) {
            // This takes a while in Spark
            numRows = (int) r.numRows();
        } else {
            numRows = cachedNumRows;
        }

        final int numCols = (int) r.numCols();

        // This cast is required, even though it would not seem necessary, at first.  Exact reason why is unknown.
        //   Will fail compilation if the cast is removed.
        final Vector [] rowVectors = (Vector []) r.rows().collect();

        final RealMatrix result = new Array2DRowRealMatrix(numRows, numCols);
        for (int i = 0; i < numRows; i++) {
            result.setRow(i, rowVectors[i].toArray() );
        }
        return result;
    }

    private static double invertSVWithTolerance(final double sv, final double tol){
        if (sv <= tol) {
            return 0;
        } else {
            return 1/sv;
        }
    }
}
