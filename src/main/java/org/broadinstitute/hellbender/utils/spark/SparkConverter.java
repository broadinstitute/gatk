package org.broadinstitute.hellbender.utils.spark;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.LinkedList;

/**
 * Class with helper methods to convert objects (mostly matrices) to/from Spark (particularly, in MLLib)
 */
public class SparkConverter {
    private static final Logger logger = LogManager.getLogger(SparkConverter.class);

    private SparkConverter() {
    }

    /**
     * Create a distributed matrix given an Apache Commons RealMatrix.
     *
     * @param sc Never {@code null}
     * @param realMat Apache Commons RealMatrix.  Never {@code null}
     * @return A distributed Spark matrix
     */
    public static RowMatrix convertRealMatrixToSparkRowMatrix(JavaSparkContext sc, RealMatrix realMat, int numSlices) {
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
        final JavaRDD<Vector> rows = sc.parallelize(rowsList, numSlices);

        // Create a RowMatrix from JavaRDD<Vector>.
        final RowMatrix mat = new RowMatrix(rows.rdd());
        logger.info("Done converting matrix to distributed Spark matrix...");
        return mat;
    }

    /**
     * Convert a local (not distributed) Spark Matrix to an Apache Commons matrix.
     *
     * @param r Never {@code null}
     * @return Not {@code null}
     */
    public static RealMatrix convertSparkMatrixToRealMatrix(final Matrix r){
        final RealMatrix result = new Array2DRowRealMatrix(r.numRows(), r.numCols());
        final double [] columnMajorMat = r.toArray();
        for (int i = 0; i < r.numRows(); i++) {
            result.setRow(i, Arrays.copyOfRange(columnMajorMat, i * r.numCols(), i * r.numCols() + r.numCols()) );
        }
        return result;
    }

    /**
     * Create an Apache Commons RealMatrix from a Spark RowMatrix.
     *
     * @param r Never {@code null}
     * @param cachedNumRows Checking the number of rows in {@code r} can be time-consuming.  Provide the value here, if it is already known.
     *                      Use {@code -1} if unknown.
     * @return Never {@code null}
     */
    public static RealMatrix convertSparkRowMatrixToRealMatrix(final RowMatrix r, final int cachedNumRows) {

        Utils.nonNull(r, "Input row matrix cannot be null");

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
}
