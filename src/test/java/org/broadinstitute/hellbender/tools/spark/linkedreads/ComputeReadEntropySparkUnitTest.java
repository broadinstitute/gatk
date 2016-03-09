package org.broadinstitute.hellbender.tools.spark.linkedreads;

import org.testng.annotations.Test;

public class ComputeReadEntropySparkUnitTest {

    @Test
    public void testComputeEntropy() {
        System.out.println(ComputeReadEntropySpark.computeEntropy("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        System.out.println(ComputeReadEntropySpark.computeEntropy("TTTTTACTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        System.out.println(ComputeReadEntropySpark.computeEntropy("ACTGCATATTCTCACTTATAAGTGGGAGCTAAATGATAAGAACTTATGAACACAAAGAAGGAAACATCAGACACTGGGGTCTAGTTGAGTGGGAAAGGTGGGAGGAGGGAGAGGAGCAGAAAAGATAATTATTGGGTACTGGGTTTAATTC"));
    }
}