package org.broadinstitute.hellbender.utils.svd;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.spark.api.java.JavaSparkContext;

/**
 *  Perform singular value decomposition (and pseudoinverse calculation).
 */
public interface SingularValueDecomposer {
    public SVD createSVD(final RealMatrix m);

    /**
     * Create the default non-spark decomposer.
     */
    public static SingularValueDecomposer getDefault(){
        return getDefault(null);
    }

    /**
     * Create the default decomposer (may be spark).
     */
    public static SingularValueDecomposer getDefault(final JavaSparkContext ctx){
        if (ctx == null) {
            return new OjAlgoSingularValueDecomposer();
        }
        return new SparkSingularValueDecomposer(ctx);
    }
}
