package org.broadinstitute.hellbender.tools.dragstr;

import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * Represents a decimation table.
 * <p>
 * A decimation table (DT) has one entry for each possible period and repeat length.
 * Each entry's value indicates the lowest "mask" bit (0-based) that is allowed to be set
 * for sampled sites.
 * </p>
 * <p>
 *     The "mask" of a site is roughly the number of sites with same period and repeat-length that
 *     have been noted before this one; actually some number in the series are skipped at the beginning
 *     of each contig (as many as the contig 0-based index in the dictionary).
 * </p>
 * <p>
 *     Thus if the DT entry for a period and repeat-length combination is 0, all sites are sampled, if it is 1
 *     every second site is discarded, if it 10, one every 1024 sites is sampled and the rest discarded.
 * </p>
 */
public class STRDecimationTable {

    /**
     * Default decimation matrix. Missing entries are assumed to be 0.
     */
    private static final int[][] DEFAULT_DECIMATION_MATRIX = new int[][] {
            {0}, // 0, 0, 0, 0, 0, 0, 0, 0 ...
            {0, 10, 10, 9, 8, 7, 5, 3, 1, 0},
            {0, 0, 9, 6, 3, 0}, // 0, 0, 0 ...
            {0, 0, 8, 4, 1, 0},
            {0, 0, 6, 0},
            {0, 0, 5, 0},
            {0, 0, 4, 0},
            {0, 0, 1, 0},
            {0}};

    /**
     * String that represents the special DT with no decimation ... ie all entries set to 0 and every site is kept.
     */
    public static final String NO_DECIMATION_STR = "NONE";

    /**
     * Strings that represents the default decimation table.
     */
    public static final String DEFAULT_DECIMATION_STR = "DEFAULT";

    /**
     * The default decimation table.
     */
    public static final STRDecimationTable DEFAULT = new STRDecimationTable(DEFAULT_DECIMATION_STR);

    /**
     * The no-decimating table.
     */
    public static final STRDecimationTable NONE = new STRDecimationTable(NO_DECIMATION_STR);

    private final String description;

    private final int[][] decimationMatrix;
    private final long[][] decimationMask;

    private final long[][] counts;

    /**
     * Creates a decimation table from its string representation. It might be a special table (NONE, DEFAULT) or the
     * path to a file that contains the table.
     */
    public STRDecimationTable(final String spec) {
        Utils.nonNull(spec);
        if (spec.equalsIgnoreCase(NO_DECIMATION_STR)) {
            decimationMatrix = new int[][] {{0}};
        } else if (spec.equalsIgnoreCase(DEFAULT_DECIMATION_STR)) {
            decimationMatrix = DEFAULT_DECIMATION_MATRIX;
        } else {
            decimationMatrix = parseDecimationMatrixFromPath(spec);
        }
        description = spec;
        decimationMask = calculateDecimationMask(decimationMatrix);
        counts = composeDecimationCounts(decimationMask);
    }

    private long[][] composeDecimationCounts(final long[][] decimationMask) {
        final long[][] result = new long[decimationMask.length][];
        for (int i = 0; i < result.length; i++) {
            result[i] = new long[decimationMask[i].length];
        }
        return result;
    }

    private static int[][] parseDecimationMatrixFromPath(String spec) {
        try (final BufferedReader reader = new BufferedReader(IOUtils.makeReaderMaybeGzipped(Paths.get(spec)))) {
            final String[][] values = reader.lines()
                    .filter(str -> !str.startsWith("#") && !str.trim().isEmpty())
                    .map(str -> Arrays.stream(str.split("\\s+"))
                               //.mapToDouble(Double::parseDouble)
                               .toArray(String[]::new))
                    .toArray(String[][]::new);
            return coerceStringMatrix(values, spec);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(spec, ex);
        } catch (final NumberFormatException ex){
            throw new UserException.BadInput(String.format("input decimation file %s contains non-valid values: %s", spec, ex.getMessage()));
        }
    }

    public void print(final PrintWriter writer) {
        Utils.nonNull(writer);
        for (final int[] row : decimationMatrix) {
            writer.println(Utils.join("\t", row));
        }
        writer.flush();
    }

    private static int[][] coerceStringMatrix(final String[][] values, final String path) {
        Utils.nonNull(values);
        if (values.length == 0) {
            LogManager.getLogger(STRDecimationTable.class)
                    .warn("Decimation matrix path provided does not seem to contain any values, we will proceed without any decimation");
            return new int[0][];
        } else {
            int totalValues = 0;
            final int[][] result = new int[values.length][];
            for (int i = 0; i < values.length; i++) {
                final String[] row = values[i];
                final int[] rowValues = new int[values[i].length];
                for (int j = 0; j < row.length; j++) {
                    final String str = row[j];
                    final int value;
                    try {
                        value = Integer.parseInt(str);
                    } catch (final NumberFormatException ex) {
                        throw badDecimationValueException(str, path, i, j, "not a valid double literal");
                    }
                    if (value < 0) {
                        throw badDecimationValueException(str, path, i, j, "negatives are not allowed");
                    } else if (Double.isNaN(value)) {
                        throw badDecimationValueException(str, path, i, j, "NaN are not allowed");
                    } else if (!Double.isFinite(value)) {
                        throw badDecimationValueException(str, path, i, j, "must be finite");
                    }
                    rowValues[j] = value;
                    totalValues++;
                }
                result[i] = rowValues;
            }
            if (totalValues == 0) {
                throw new UserException.BadInput("the input decimation matrix does contain any values:" + path);
            }
            return result;
        }
    }

    private static RuntimeException badDecimationValueException(final String str, final String path, final int i, final int j,
                                                                final String details) {
        throw new UserException.BadInput(String.format("bad decimation value found in %s for period and repeats (%d, %d) with string (%s)%s",
                path, i, j, str, details == null || details.isEmpty()? "": ": " + details));
    }

    private static long[][] calculateDecimationMask(final int[][] decimationMatrix) {
        Utils.nonNull(decimationMatrix);
        final long[][] result = new long[decimationMatrix.length][];
        for (int i = 0; i < result.length; i++) {
            final int[] row = decimationMatrix[i];
            result[i] = new long[row.length];
            for (int j = 0; j < row.length; j++) {
                result[i][j] = (1 << row[j]) - 1;
            }
        }
        return result;
    }

    public long mask(final int period, final int repeats) {
        final int p = period >= decimationMask.length ? decimationMask.length - 1 : period;
        final long[] masks = decimationMask[p];
        if (masks.length == 0) {
            return 0;
        } else if (repeats >= masks.length) {
            return masks[masks.length - 1];
        } else {
            return masks[repeats];
        }
    }

    public boolean decimate(final long mask, final int bestPeriod, final long bestPeriodRepeats) {
        if (counts.length <= bestPeriod) {
            return false;
        } else {
            final long[] periodCounts = counts[bestPeriod];
            if (bestPeriodRepeats >= periodCounts.length) {
                return false;
            } else {
                final long left = mask;
                final long right = decimationMask[bestPeriod][(int) bestPeriodRepeats];
                return ((int) left & (int) right) != 0 || ((left >> 32) & (right >> 32)) != 0;
            }
        }
    }

    public int decimationBit(int period, int repeatCount) {
        if (period >= decimationMatrix.length) {
            return 0;
        } else if (repeatCount >= decimationMatrix[period].length) {
            return 0;
        } else {
            return decimationMatrix[period][repeatCount];
        }
    }

    @Override
    public String toString() {
        return description;
    }
}
