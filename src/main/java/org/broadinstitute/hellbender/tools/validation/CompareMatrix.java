package org.broadinstitute.hellbender.tools.validation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;

/**
 * CompareMatrix contains a square matrix of linear dimension QualityUtils.MAX_SAM_QUAL_SCORE.
 * The elements are counts of the number of times bam 1 had a quality score of i and bam 2 had a quality
 * score of j, (stored in mat[i][j]).
 */
public final class CompareMatrix implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final int dimension = QualityUtils.MAX_SAM_QUAL_SCORE+1;
    private final byte[] bin;

    private long[][] mat;
    private long[][] binnedMat;

    public CompareMatrix(final byte[] binning) {
        Utils.nonNull(binning);
        this.bin = binning.clone();
        mat = new long[dimension][dimension];
        binnedMat = new long[dimension][dimension];
    }

    public void add(final byte[] first, final byte[] second) {
        if (first.length != second.length) {
            throw new UserException("The length of the quality scores are not the same for read " + first.length + "," + second.length);
        }

        for (int i = 0; i < first.length; ++i) {
            mat[first[i]][second[i]]++;
            binnedMat[bin[first[i]]][bin[second[i]]]++;
        }
    }

    public CompareMatrix add(final CompareMatrix c) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                mat[i][j] += c.mat[i][j];
                binnedMat[i][j] += c.binnedMat[i][j];//only need to bin when adding reads, not here
            }
        }
        return this;
    }

    public boolean hasNonDiagonalElements() {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (i != j && mat[i][j] != 0) {
                    return true;
                }
            }
        }
        return false;
    }

    private static void printSummary(final PrintStream printStream, final long[][] matrix, final String name) {
        final int totalSize = 2*(dimension-1)+1;
        final long[] deltaColumns = new long[totalSize];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                deltaColumns[(dimension-1) + i-j] += matrix[i][j];
            }
        }
        printStream.println("-----------" + name + " summary------------");
        final long columnsSum = MathUtils.sum(deltaColumns);
        if (columnsSum == deltaColumns[dimension - 1]) {
            printStream.println("all " + columnsSum + " quality scores are the same");
        } else {
            printStream.println("diff" + "\t" + "count" + "\t" + "%total");
            for (int k = 0; k < totalSize; ++k) {
                if (deltaColumns[k] != 0) {
                    printStream.println(String.format("%d\t%d\t%.4f", k - (dimension - 1), deltaColumns[k], deltaColumns[k]*100.0/columnsSum));
                }
            }
        }

    }

    /**
     *  Prints the matrix of counts of quality score from read 1 vs quality score from read 2.
     *  Only prints the non-zero entries for readability.
     */
    private static void print(final PrintStream printStream, final long[][] matrix, final String name) {
        printStream.println("---------" + name + " full matrix (non-zero entries) ----------");

        printStream.println("QRead1" + "\t" + "QRead2" + "\t" + "diff" + "\t" + "count");
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (matrix[i][j] != 0){
                    printStream.println(i + "\t" + j + "\t" + (i-j) + "\t" + matrix[i][j]);
                }
            }
        }
    }

    public void printOutResults(final String outputFilename) {
        if (outputFilename != null) {
            try (final OutputStream os = new FileOutputStream(outputFilename)) {
                final PrintStream ps = new PrintStream(os);
                this.printOutput(ps);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile("unable to write to output file: " + outputFilename, e);
            }
        } else {
            this.printOutput(System.out);
        }
    }

    private void printOutput(final PrintStream ps) {
        printSummary(ps, mat, "CompareMatrix");
        ps.println();
        print(ps, mat, "CompareMatrix");
        printSummary(ps, binnedMat, "CompareMatrix-binned");
        ps.println();
        print(ps, binnedMat, "CompareMatrix-binned");
    }
}
