package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntegerCopyNumberTransitionMatrixCollectionUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/germlinehmm";
    private static final File HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE = new File(TEST_SUB_DIR,
            "homo_sapiens_germline_HMM_priors.tsv");
    private final Set<String> HOMO_SAPIENS_SEX_GENOTYPES = Arrays.stream(new String[] {"SEX_XX", "SEX_XY"})
            .collect(Collectors.toSet());
    private final Set<String> HOMO_SAPIENS_ALL_CONTIGS = Arrays.stream(new String[] {"1", "2", "3", "4", "5", "6",
        "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"})
            .collect(Collectors.toSet());
    private final Set<String> HOMO_SAPIENS_AUTOSOMAL_CONTIGS = Arrays.stream(new String[] {"1", "2", "3", "4", "5", "6",
            "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"})
            .collect(Collectors.toSet());

    @Test
    public void testSucessfulLoading() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        collection.assertCompleteness(HOMO_SAPIENS_SEX_GENOTYPES, HOMO_SAPIENS_ALL_CONTIGS);
    }

    @Test
    public void testValidity() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        /* autosomal contigs */
        for (final String sexGenotype : HOMO_SAPIENS_SEX_GENOTYPES) {
            for (final String contig : HOMO_SAPIENS_AUTOSOMAL_CONTIGS) {
                IntegerCopyNumberTransitionMatrixDataUnitTest.assertRealMatrixEquals(
                        IntegerCopyNumberTransitionMatrixDataUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_AUTOSOMAL_TRUTH,
                        collection.get(sexGenotype, contig).getTransitionMatrix(), 1e-16);
            }
        }
        /* allosomal contigs */
        IntegerCopyNumberTransitionMatrixDataUnitTest.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixDataUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_X_TRUTH,
                collection.get("SEX_XX", "X").getTransitionMatrix(), 1e-16);
        IntegerCopyNumberTransitionMatrixDataUnitTest.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixDataUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_Y_TRUTH,
                collection.get("SEX_XX", "Y").getTransitionMatrix(), 1e-16);
        IntegerCopyNumberTransitionMatrixDataUnitTest.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixDataUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_X_TRUTH,
                collection.get("SEX_XY", "X").getTransitionMatrix(), 1e-16);
        IntegerCopyNumberTransitionMatrixDataUnitTest.assertRealMatrixEquals(
                IntegerCopyNumberTransitionMatrixDataUnitTest.HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_Y_TRUTH,
                collection.get("SEX_XY", "Y").getTransitionMatrix(), 1e-16);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testMissingSexGenotypeAssertionFailure() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        collection.get("SEX_XYZ", "1");
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testMissingContigAssertionFailure() {
        final IntegerCopyNumberTransitionMatrixCollection collection =
                IntegerCopyNumberTransitionMatrixCollection.read(HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE);
        collection.get("SEX_XY", "24");
    }
}
