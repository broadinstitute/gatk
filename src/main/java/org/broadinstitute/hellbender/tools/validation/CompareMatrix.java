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

    private long[][] mat;

    public CompareMatrix() {
        mat = new long[dimension][dimension];
    }

    public void add(byte[] first, byte[] second) {
        if (first.length != second.length) {
            throw new UserException("The length of the quality scores are not the same for read " + first.length + "," + second.length);
        }

        for (int i = 0; i < first.length; ++i) {
            mat[(int) first[i]][(int) second[i]]++;
        }
    }

    public CompareMatrix add(CompareMatrix c) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                mat[i][j] += c.mat[i][j];
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

    public void printSummary(final PrintStream printStream) {
        Utils.nonNull(printStream);
        final int totalSize = 2*(dimension-1)+1;
        final long[] deltaColumns = new long[totalSize];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                deltaColumns[(dimension-1) + i-j] += mat[i][j];
            }
        }
        printStream.println("-----------CompareMatrix summary------------");
        if (MathUtils.sum(deltaColumns) == deltaColumns[dimension - 1]) {
            printStream.println("all " + MathUtils.sum(deltaColumns) + " quality scores are the same");
        } else {
            printStream.println("diff" + "\t" + "count");
            for (int k = 0; k < totalSize; ++k) {
                if (deltaColumns[k] != 0) {
                    printStream.println(k - (dimension - 1) + "\t" + deltaColumns[k]);
                }
            }
        }

    }

    /**
     *  Prints the matrix of counts of quality score from read 1 vs quality score from read 2.
     *  Only prints the non-zero entries for readability.
     */
    public void print(final PrintStream printStream) {
        Utils.nonNull(printStream);
        printStream.println("---------CompareMatrix full matrix (non-zero entries) ----------");

        printStream.println("QRead1" + "\t" + "QRead2" + "\t" + "diff" + "\t" + "count");
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (mat[i][j] != 0){
                    printStream.println(i + "\t" + j + "\t" + (i-j) + "\t" + mat[i][j]);
                }
            }
        }
    }


    public void printOutput(final String outputFilename) {
        if (outputFilename != null) {
            try (OutputStream os = new FileOutputStream(outputFilename)) {
                final PrintStream ps = new PrintStream(os);
                printSummary(ps);
                ps.println();
                print(ps);
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile("unable to write to output file: " + outputFilename, e);
            }
        } else {
            // Print to stdout instead
            final PrintStream ps = System.out;
            printSummary(ps);
            ps.println();
            print(ps);
        }
    }
}
