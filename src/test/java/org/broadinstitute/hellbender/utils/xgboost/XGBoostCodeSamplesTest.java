package org.broadinstitute.hellbender.utils.xgboost;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;
import ml.dmlc.xgboost4j.java.DMatrix;

public class XGBoostCodeSamplesTest extends GATKBaseTest {

    @Test
    public void testEmpty() {

    }

    // docs on https://xgboost.readthedocs.io/en/stable/jvm/java_intro.html
    @Test
    public  void testOne() throws Exception {
        float[] data = new float[] {1f,2f,3f,4f,5f,6f};
        int nrow = 3;
        int ncol = 2;
        float missing = 0.0f;
        DMatrix dmat = new DMatrix(data, nrow, ncol, missing);
        float[] weights = new float[] {1f,2f,1f};
        dmat.setWeight(weights);

    }




}
