package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * Entry point for creating an instance of SVD.  When the object is created, all of the calculation will be done as well.
 */
public final class SVDFactory {

    /**
     * Create a SVD instance using Apache Commons Math.
     *
     * @param m matrix that is not {@code null}
     * @return SVD instance that is never {@code null}
     */
    public static SVD createSVD(final RealMatrix m){
        Utils.nonNull(m, "Cannot perform SVD on a null.");
        return createSVD(m, null);
    }

    /**
     * Create a SVD instance using a spark context.
     *
     * @param m matrix that is not {@code null}
     * @param ctx JavaSparkContext.  {@code null} is allowed, but will fall back to Apache Commons Math implementation.
     * @return SVD instance that is never {@code null}
     */
    public static SVD createSVD(final RealMatrix m, final JavaSparkContext ctx){
        Utils.nonNull(m, "Cannot create SVD from a null matrix.");
        if (ctx == null) {
            return new OjAlgoSingularValueDecomposer().createSVD(m);
        }
        return new SparkSingularValueDecomposer(ctx).createSVD(m);
    }
}
