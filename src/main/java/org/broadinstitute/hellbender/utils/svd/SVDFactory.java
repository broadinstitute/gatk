package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Entry point for creating an instance of SVD.  When the object is created, all of the calculation will be done as well.
 *
 */
public class SVDFactory {

    /** Create a SVD instance using Apache Commons Math.
     *
     * @param m matrix that is not {@code null}
     * @return SVD instance that is never {@code null}
     */
    public static SVD createSVD(final RealMatrix m){
        Utils.nonNull(m, "Cannot perform SVD on a null.");
        return createSVD(m, null);
    }

    /**
     * Create a SVD instance using A spark context.
     *
     * @param m matrix that is not {@code null}
     * @param ctx JavaSparkContext.  {@code null} is allowed, but will fall back to Apache Commons Math implementation.
     * @return SVD instance that is never {@code null}
     */
    public static SVD createSVD(final RealMatrix m, final JavaSparkContext ctx){
        if (ctx == null) {
            return ApacheSingularValueDecomposer.createSVD(m);
        }
        return SparkSingularValueDecomposer.createSVD(ctx, m);
    }
}
