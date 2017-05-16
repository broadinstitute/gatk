package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.utils.MathObjectAsserts;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntegerCopyNumberTransitionMatrixCollectionUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/germlinehmm";
    private static final File HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE = new File(TEST_SUB_DIR,
            "homo_sapiens_germline_HMM_priors.tsv");

    @Test
    public void testSuccessfulLoading() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        collection.assertCompleteness();
    }

    @Test
    public void testValidity() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        /* autosomal contigs */
        for (final String sexGenotype : IntegerCopyNumberTransitionProbabilityCacheCollectionUnitTest.HOMO_SAPIENS_SEX_GENOTYPES) {
            for (final String contig : IntegerCopyNumberTransitionProbabilityCacheCollectionUnitTest.HOMO_SAPIENS_AUTOSOMAL_CONTIGS) {
                MathObjectAsserts.assertRealMatrixEquals(
                        IntegerCopyNumberTransitionMatrixUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_AUTOSOMAL_TRUTH,
                        collection.get(sexGenotype, contig).getTransitionMatrix());
            }
        }
        /* allosomal contigs */
        MathObjectAsserts.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_X_TRUTH,
                collection.get("SEX_XX", "X").getTransitionMatrix());
        MathObjectAsserts.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_Y_TRUTH,
                collection.get("SEX_XX", "Y").getTransitionMatrix());
        MathObjectAsserts.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_X_TRUTH,
                collection.get("SEX_XY", "X").getTransitionMatrix());
        MathObjectAsserts.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_Y_TRUTH,
                collection.get("SEX_XY", "Y").getTransitionMatrix());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingSexGenotypeAssertionFailure() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        collection.get("SEX_XYZ", "1");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMissingContigAssertionFailure() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        collection.get("SEX_XY", "24");
    }
}
