package org.broadinstitute.hellbender.tools.spark.validation;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Serializable;

/**
 * CompareMatrix contains a square matrix of linear dimension QualityUtils.MAX_SAM_QUAL_SCORE.
 * The elements are counts of the number of times bam 1 had a quality score of i and bam 2 had a quality
 * score of j, (stored in mat[i][j]).
 */
class CompareMatrix implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final int dimension = QualityUtils.MAX_SAM_QUAL_SCORE+1;

    private int[][] mat;

    CompareMatrix() {
        mat = new int[dimension][dimension];
    }

    void add(byte[] first, byte[] second) {
        if (first.length != second.length) {
            throw new UserException("The length of the quality scores are not the same for read " + first.length + "," + second.length);
        }

        for (int i = 0; i < first.length; ++i) {
            mat[(int) first[i]][(int) second[i]]++;
        }
    }

    CompareMatrix add(CompareMatrix c) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                mat[i][j] += c.mat[i][j];
            }
        }
        return this;
    }

    Boolean hasNonDiagonalElements() {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if (i != j && mat[i][j] != 0) {
                    return true;
                }
            }
        }
        return false;
    }

    void printSummary(OutputStream outputStream) {
        PrintStream printStream = new PrintStream(outputStream);
        // This prints a two-column summary with the difference in quality scores between the two BAMs.
        // If there was only one read per bam with [10,11,12,13] and [9,11,13,13], the summary would be
        // |-93|-92|...|-1|0|1|...|92|93|
        // |0  |0  |...|1 |2|1|...|0 |0 |

        printStream.println("-----------CompareMatrix summary------------");
        int totalSize = 2*(dimension-1)+1;
        int[] deltaColumns = new int[totalSize];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                deltaColumns[(dimension-1) + i-j] += mat[i][j];
            }
        }
        int[] maxColumnWidths = new int[totalSize];
        for (int k = 0; k < totalSize; ++k) {
            int size =  Integer.toString(k - (dimension-1)).length();
            maxColumnWidths[k] = size;
            // Increase the column width if the with of the count is larger than the label width.
            size = Integer.toString(deltaColumns[k]).length();
            if (size > maxColumnWidths[k]) {
                maxColumnWidths[k] = size;
            }
        }
        for (int k = 0; k < totalSize; ++k) {
            printStream.print("|");
            String fmt = "%-" + Integer.toString(maxColumnWidths[k]) + "s";
            printStream.format(fmt, k - (dimension-1));
        }
        printStream.println("|");

        for (int k = 0; k < totalSize; ++k) {
            printStream.print("|");
            String fmt = "%-" + Integer.toString(maxColumnWidths[k]) + "s";
            printStream.format(fmt, deltaColumns[k]);
        }
        printStream.println("|");
    }

    void print(OutputStream outputStream) {
        PrintStream printStream = new PrintStream(outputStream);
        // Prints the full matrix of counts of quality score from read 1 vs quality score from read 2.
        printStream.println("---------CompareMatrix full matrix----------");

        // First find the width (in chars) of each column (for pretty printing)
        int[] maxColumnWidths = new int[dimension];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                int size = Integer.toString(mat[i][j]).length();
                if (size > maxColumnWidths[j]) {
                    maxColumnWidths[j] = size;
                }
            }
        }

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                printStream.print("|");
                String fmt = "%-" + Integer.toString(maxColumnWidths[j]) + "s";
                printStream.format(fmt, mat[i][j]);
                if (j == dimension - 1) {
                    printStream.println("|");
                }
            }
        }
    }
}
