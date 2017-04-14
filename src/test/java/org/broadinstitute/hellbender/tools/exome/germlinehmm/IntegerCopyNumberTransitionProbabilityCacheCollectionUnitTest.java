package org.broadinstitute.hellbender.tools.exome.germlinehmm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class IntegerCopyNumberTransitionProbabilityCacheCollectionUnitTest extends BaseTest {
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

    /**
     * The maximum copy number for all transition matrices in this test suite
     */
    private final int MAX_COPY_NUMBER = 10;

    private final int[] DISTANCES = {1, 2, 5, 10, 100, 500, 1000};

    private static final double EPSILON = 1e-10;

    @Test
    public void performCompleteTestNoPadding() {
        final IntegerCopyNumberTransitionProbabilityCacheCollection cache =
                new IntegerCopyNumberTransitionProbabilityCacheCollection(
                        HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE, false);
        performCompleteTest(cache, false);
    }

    @Test
    public void performCompleteTestWithPadding() {
        final IntegerCopyNumberTransitionProbabilityCacheCollection cache =
                new IntegerCopyNumberTransitionProbabilityCacheCollection(
                        HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE, true);
        performCompleteTest(cache, true);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMaxCopyNumberNoPadding() {
        final IntegerCopyNumberTransitionProbabilityCacheCollection cache =
                new IntegerCopyNumberTransitionProbabilityCacheCollection(
                        HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE, false);
        cache.getMaxCopyNumber();
    }

    @Test
    public void testMaxCopyNumberWithPadding() {
        final IntegerCopyNumberTransitionProbabilityCacheCollection cache =
                new IntegerCopyNumberTransitionProbabilityCacheCollection(
                        HOMO_SAPIENS_COPY_NUMBER_TRANSITION_PRIOR_TABLE_FILE, true);
        Assert.assertEquals(cache.getMaxCopyNumber(), MAX_COPY_NUMBER);
    }

    private void performCompleteTest(final IntegerCopyNumberTransitionProbabilityCacheCollection cache,
                                     final boolean padded) {

        for (final String sexGenotype : HOMO_SAPIENS_SEX_GENOTYPES) {
            for (final int dist : DISTANCES) {
                for (final String contig : HOMO_SAPIENS_ALL_CONTIGS) {

                    /* set the per-base transition matrix according to contig and sex genotype */
                    final RealMatrix perBaseTransitionMatrix;
                    if (HOMO_SAPIENS_AUTOSOMAL_CONTIGS.contains(contig)) {
                        perBaseTransitionMatrix = IntegerCopyNumberTransitionMatrixDataUnitTest
                                .HOMO_SAPIENS_COPY_NUMBER_TRANSITION_AUTOSOMAL_TRUTH;
                    } else if (contig.equals("X")) {
                        if (sexGenotype.equals("SEX_XX")) {
                            perBaseTransitionMatrix = IntegerCopyNumberTransitionMatrixDataUnitTest
                                    .HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_X_TRUTH;
                        } else { /* SEX_XY */
                            perBaseTransitionMatrix = IntegerCopyNumberTransitionMatrixDataUnitTest
                                    .HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_X_TRUTH;
                        }
                    } else { /* contig = Y */
                        if (sexGenotype.equals("SEX_XX")) {
                            perBaseTransitionMatrix = IntegerCopyNumberTransitionMatrixDataUnitTest
                                    .HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XX_Y_TRUTH;
                        } else { /* SEX_XY */
                            perBaseTransitionMatrix = IntegerCopyNumberTransitionMatrixDataUnitTest
                                    .HOMO_SAPIENS_COPY_NUMBER_TRANSITION_XY_Y_TRUTH;
                        }
                    }

                    /* the list of copy number states */
                    final List<IntegerCopyNumberState> integerCopyNumberStates =
                            IntStream.range(0, cache.getMaxCopyNumber(sexGenotype, contig) + 1)
                                    .mapToObj(IntegerCopyNumberState::new)
                                    .collect(Collectors.toList());

                    /* calculate the log transition matrix directly */
                    final RealMatrix transitionMatrixDirect = getDirectMatrixPower(perBaseTransitionMatrix, dist);

                    /* calculate the log transition matrix using the cache class */
                    final RealMatrix transitionMatrixFromCache = getTransitionMatrix(dist, sexGenotype, contig,
                            integerCopyNumberStates, integerCopyNumberStates,  cache);

                    /* if the collection is padded, we must pad the truth as well */
                    final RealMatrix transitionMatrixDirectPadded;
                    if (padded) {
                        transitionMatrixDirectPadded = MatrixUtils.createRealIdentityMatrix(MAX_COPY_NUMBER + 1);
                        for (int i = 0; i < transitionMatrixDirect.getRowDimension(); i++) {
                            for (int j = 0; j < transitionMatrixDirect.getColumnDimension(); j++) {
                                transitionMatrixDirectPadded.setEntry(i, j, transitionMatrixDirect.getEntry(i, j));
                            }
                        }
                    } else {
                        transitionMatrixDirectPadded = transitionMatrixDirect;
                    }
                    assertEqualMatrices(transitionMatrixDirectPadded, transitionMatrixFromCache);
                }
            }
        }
    }

    private static RealMatrix getDirectMatrixPower(final RealMatrix mat, final int power) {
        if (power < 0) {
            throw new IllegalArgumentException("Can not calculate negative matrix powers");
        } else if (power == 0) {
            return MatrixUtils.createRealIdentityMatrix(mat.getColumnDimension());
        } else {
            return mat.power(power);
        }
    }

    private static RealMatrix getTransitionMatrix(final int distance,
                                                  final String sexGenotype, final String contig,
                                                  final List<IntegerCopyNumberState> toList,
                                                  final List<IntegerCopyNumberState> fromList,
                                                  final IntegerCopyNumberTransitionProbabilityCacheCollection cache) {
        final RealMatrix mat = new Array2DRowRealMatrix(toList.size(), fromList.size());
        for (int i = 0; i < toList.size(); i++) {
            for (int j = 0; j < fromList.size(); j++) {
                mat.setEntry(i, j, FastMath.exp(cache.logTransitionProbability(distance, sexGenotype, contig,
                        toList.get(i), fromList.get(j))));
            }
        }
        return mat;
    }

    private static void assertEqualMatrices(final RealMatrix mat1, final RealMatrix mat2) {
        Assert.assertEquals(mat1.subtract(mat2).getNorm(), 0, EPSILON);
    }
}
