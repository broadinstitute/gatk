package org.broadinstitute.hellbender.utils.dragstr;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.dragstr.DragstrHyperParameters;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Holds the values of the DRAGstr model parameters for different combinations of repeat unit length (period) and
 * number of repeat units.
 * <p>
 *     This class is immutable.
 * </p>
 * <h3>Parameter matrices</h3>
 * <p>
 *     The DRAGstr model are composed of three matrices with the same dimensions (period x repeat-length):
 *     <dl>
 *         <dt>GOP</dt>
 *         <dd>gap-opening-penalty table where we store the Phred scale score that we use for pairhmm match to insert or
 *         delete transitions</dd>.
 *         <dt>GCP</dt>
 *         <dd>gap-continuation-penalty table where we store the Phred scale score that we use for pair-hmm insert to insert, and
 *             deletion to deletion
 *             </dd>.
 *         <dt>API</dt>
 *         <dd>adjusted prior probability of an indel occurring at a position in the reference with an STR.</dd>
 *     </dl>
 * </p>
 * <h3>Text file format</h3>
 * If P is the maximum period length and L is the maximum repeat-length:
 * <p>
 *     <pre>
 *         ############################
 *         # DragstrParams
 *         # -------------------------
 *         # annotation1 = value1
 *         # annotation2 = value2
 *         # ...
 *         ###########################
 *              1      2     3     4 ...     L
 *         GOP:
 *          gop11  gop12 gop13 gop14 ... gop1L
 *          gop21  gop22 ...
 *          ...
 *          gopP1  gopP2 ...             gopPL
 *         GCP:
 *          gcp11  gcp12 gcp13 ...
 *          ...
 *          gcpP1  apiP2 ...     ...     apiPL
 *         API:
 *          api11  api12 ...
 *          ..
 *          apiP1  apiP2 ...             apiPL
 *     </pre>
 * </p>
 * <p>
 *     Currently there is not a fix/standard set of annotations but is just a mechanism to add some additional
 *     metadata that might be useful.
 * </p>
 */
public final class DragstrParams {

    public static final String GOP_TABLE_NAME = "GOP";
    public static final String GCP_TABLE_NAME = "GCP";
    public static final String API_TABLE_NAME = "API";

    /**
     * String returned by {@link #toString()} when this object has not been retrieved from nor persisted into a file.
     */
    private static final String NOT_PERSISTENT_STRING = "<non-persistent>";

    /**
     * Initial size of the buffer/builder used to compose lines in the output text format.
     */
    private static int LINE_BUILDER_BUFFER_SIZE = 1024;

    /**
     * Holds the path of the file this param where loaded from or persistent into last.
     */
    private String path;

    /**
     * Gap-Open-Penalty default Phred scores.
     * <p>
     *     Each row represent a different period (index 0 is period length == 1, index 1 is period length == 2) to the default maximum of {@value DragstrHyperParameters#DEFAULT_MAX_PERIOD}.
     * </p>
     * <p>
     *     Then each element (column) in that row is the value for the ith repeat length in units (0-based) up to the default maximum of {@value DragstrHyperParameters#DEFAULT_MAX_REPEAT_LENGTH}.
     * </p>
     *
     * <p>
     *     This values were copied from DRAGstr documentation provided by Illumina.
     *     It is unknown to us how they came out with this values nor we know of any way to generate them
     *     on the fly.
     * </p>
     */
    //@formatter:off -- this disables IntelliJ code reformatting if preferences are set; Preferences > Editor > Code Style > Formatter Control
    private static final double[][] DEFAULT_GOP = { /* GOP */
     //   Repeat Length
     //             1,     2,     3,     4,     5,     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,    20
     // Period
     /* 1 */   {45.00, 45.00, 45.00, 45.00, 45.00, 45.00, 40.50, 33.50, 28.00, 24.00, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75},
     /* 2 */   {39.50, 39.50, 39.50, 39.50, 36.00, 30.00, 27.25, 25.00, 24.25, 24.75, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.75},
     /* 3 */   {38.50, 41.00, 41.00, 41.00, 41.00, 37.50, 35.25, 34.75, 34.75, 33.25, 33.25, 33.25, 32.50, 30.75, 28.50, 29.00, 29.00, 29.00, 29.00, 29.00},
     /* 4 */   {37.50, 39.00, 39.00, 37.75, 34.00, 34.00, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 31.75, 31.75, 31.75, 31.75, 31.75},
     /* 5 */   {37.00, 40.00, 40.00, 40.00, 36.00, 35.00, 24.50, 24.50, 24.50, 24.50, 22.50, 22.50, 22.50, 23.50, 23.50, 23.50, 23.50, 23.50, 23.50, 23.50},
     /* 6 */   {36.25, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00},
     /* 7 */   {36.00, 40.50, 40.50, 40.50, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75},
     /* 8 */   {36.25, 39.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75}};
     //@formatter:on

    /**
     * API default Phred scores.
     * <p>
     *     Each row represent a different period (index 0 is period length == 1, index 1 is period length == 2) to the default maximum of {@value DragstrHyperParameters#DEFAULT_MAX_PERIOD}.
     * </p>
     * <p>
     *     Then each element (column) in that row is the value for the ith repeat length in units (0-based) up to the default maximum of {@value DragstrHyperParameters#DEFAULT_MAX_REPEAT_LENGTH}.
     * </p>
     *
     * <p>
     *     This values were copied from DRAGstr documentation provided by Illumina.
     *     It is unknown to us how they came out with this values nor we know of any way to generate them
     *     on the fly.
     * </p>
     */
    //@formatter:off -- this disables IntelliJ code reformatting if preferences are set; Preferences > Editor > Code Style > Formatter Control
    private static final double[][] DEFAULT_API = { /* API */
     //   Repeat Length
     //           1,     2,     3,     4,     5,     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,    18,    19,    20
     // Period
     /* 1 */ {39.00, 39.00, 37.00, 35.00, 32.00, 26.00, 20.00, 16.00, 12.00, 10.00,  8.00,  7.00,  7.00,  6.00,  6.00,  5.00,  5.00,  4.00,  4.00,  4.00},
     /* 2 */ {30.00, 30.00, 29.00, 22.00, 17.00, 14.00, 11.00,  8.00,  6.00,  5.00,  4.00,  4.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  2.00,  2.00},
     /* 3 */ {27.00, 27.00, 25.00, 18.00, 14.00, 12.00,  9.00,  7.00,  5.00,  4.00,  3.00,  3.00,  3.00,  3.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00},
     /* 4 */ {27.00, 27.00, 18.00,  9.00,  9.00,  9.00,  9.00,  3.00,  3.00,  3.00,  3.00,  3.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00},
     /* 5 */ {29.00, 29.00, 18.00,  8.00,  8.00,  8.00,  4.00,  3.00,  3.00,  3.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00,  2.00},
     /* 6 */ {25.00, 25.00, 10.00, 10.00, 10.00,  4.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00},
     /* 7 */ {21.00, 21.00, 11.00, 11.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00,  5.00},
     /* 8 */ {18.00, 18.00, 10.00,  6.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00}};
    //@formatter:on

    /**
     * Gap-Continuation-Penalty (GCP)
     * <p>
     *     Each row represent a different period (index 0 is period length == 1, index 1 is period length == 2) to the default maximum of {@value DragstrHyperParameters#DEFAULT_MAX_PERIOD}.
     * </p>
     * <p>
     *     Then each element (column) in that row is the value for the ith repeat length in units (0-based) up to the default maximum of {@value DragstrHyperParameters#DEFAULT_MAX_REPEAT_LENGTH}.
     * </p>
     *
     * <p>
     *    This value is not estimated and it rather has fixed values. They only depend on the period length.
     *    With a very simple formula of 10.0 / period. Since DRAGEN tables are approximated just the second
     *    decimal. We do that below using: {@code round(10 * 100 / period) / 100}
     *    It is unknown to us how they came out with this values nor we know of any way to generate them
     *    on the fly.
     * </p>
     * <p>
     *     As a result the first row (period 1, index 0) is just an array of 10.0, the second row (period 2, index 1) is filled with 5.0,
     *     the third 3.33, the fourth 2.25 and so forth.
     * </p>
     */
    private static final double[][] DEFAULT_GCP = IntStream.rangeClosed(1, DragstrHyperParameters.DEFAULT_MAX_PERIOD)
            .mapToDouble(period -> Math.round(1000.0 / period) / 100.0) // This way we end-up with two decimal precision.
            .mapToObj(value -> MathUtils.doubles(DragstrHyperParameters.DEFAULT_MAX_REPEAT_LENGTH, value))
            .toArray(double[][]::new);

    /**
     * Default parameters when there is not enough data for estimation.
     */
    public static final DragstrParams DEFAULT = new DragstrParams(
                   DragstrHyperParameters.DEFAULT_MAX_PERIOD,
                   DragstrHyperParameters.DEFAULT_MAX_REPEAT_LENGTH,
                   DEFAULT_GOP,
                   DEFAULT_GCP,
                   DEFAULT_API);

    private final int maxPeriod;
    private final int maxRepeats;
    private final double[][] gop;
    private final double[][] gcp;
    private final double[][] api;

    private Map<Object, AlleleFrequencyCalculator> afcs;

    public DragstrParams(final String path) {
        this(openBufferedReader(path), path);
    }

    private static BufferedReader openBufferedReader(String path) {
        try {
            return Files.newBufferedReader(Paths.get(path));
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(path, ex);
        }
    }

    public static DragstrParams of(final int maxPeriod, final int maxRepeats, final double[][] gop, final double[][] gcp, final double[][] api) {
        ParamUtils.isPositive(maxPeriod, "max period must be a positive");
        ParamUtils.isPositive(maxRepeats, "max repeats must be a positive");
        Utils.nonNull(gop, "gop cannot be null");
        Utils.nonNull(gcp, "gcp cannot be null");
        Utils.nonNull(api, "api cannot be null");
        Utils.validate(gop.length == maxPeriod, "input gop length must match maxPeriod");
        Utils.validate(gcp.length == maxPeriod, "input gcp length must match maxPeriod");
        Utils.validate(api.length == maxPeriod, "input api length must match maxPeriod");
        for (int i = 0; i < maxPeriod; i++) {
            final double[] gopRow = gop[i];
            final double[] gcpRow = gcp[i];
            final double[] apiRow = api[i];
            for (int j = 0; j < maxRepeats; j++) {
                Utils.validate(gopRow[j] >= 0 && Double.isFinite(gopRow[j]), "bad gop value: " + gopRow[j]);
                Utils.validate(gcpRow[j] >= 0 && Double.isFinite(gcpRow[j]), "bad gcp value: " + gcpRow[j]);
                Utils.validate(apiRow[j] >= 0 && Double.isFinite(apiRow[j]), "bad api value: " + apiRow[j]);
            }
        }
        return new DragstrParams(maxPeriod, maxRepeats, gop.clone(), gcp.clone(), api.clone());
    }

    private BufferedWriter openBufferedWriter(final String path) {
        try {
            return Files.newBufferedWriter(Paths.get(path));
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(path, ex);
        }
    }

    private DragstrParams(final BufferedReader reader, final String path) {
        try {
            this.path = path;
            String header;
            while ((header = reader.readLine()) != null) {
                if (!header.startsWith("#")) {
                    break;
                }
            }
            if (header == null) {
                throw new UserException.BadInput("there is no content in the dragstr-params file " + path);
            }
            final String[] headerParts = header.split("\\s+");
            final int[] repeats = Arrays.stream(headerParts)
                    .filter(str -> !str.isEmpty())
                    .mapToInt(str -> {
                        try {
                            return Integer.parseInt(str);
                        } catch (final NumberFormatException ex) {
                            throw new UserException.BadInput("bad format for an integer", ex);
                        }
                    })
                    .toArray();
            final int maxRepeats = repeats.length;
            for (int i = 0; i < repeats.length; i++) {
                if (repeats[i] != i + 1) {
                    throw new UserException.BadInput("the DRAGstr parameter file header line must contain integers starting at 1 " + Arrays.toString(repeats));
                }
            }
            final Map<String, double[][]> tables = new HashMap<>();
            String line = reader.readLine();
            if (line == null) {
                throw new UserException.BadInput("end of table list before expected");
            }
            String tableName = line.replaceAll(":$", "");
            List<String> tableLines = new ArrayList<>();
            while ((line = reader.readLine()) != null) {
                if (line.charAt(line.length() - 1) == ':') {
                    tables.put(tableName, linesToMatrix(tableLines, repeats.length));
                    tableName = line.replaceAll(":$", "");
                    tableLines.clear();
                } else {
                    tableLines.add(line);
                }
            }
            if (tableName.isEmpty()) {
                throw new UserException.BadInput("table with no name");
            }
            tables.put(tableName, linesToMatrix(tableLines, repeats.length));
            final double[][] gopMatrix = mandatoryMatrix(tables, GOP_TABLE_NAME);
            final double[][] gcpMatrix = mandatoryMatrix(tables, GCP_TABLE_NAME);
            final double[][] apiMatrix = mandatoryMatrix(tables, API_TABLE_NAME);
            final int maxPeriod = gopMatrix.length;
            checkMatricesAreValid(maxPeriod, maxRepeats, gopMatrix, gcpMatrix, apiMatrix);

            this.maxPeriod = maxPeriod;
            this.maxRepeats = maxRepeats;
            this.gop = gopMatrix;
            this.gcp = gcpMatrix;
            this.api = apiMatrix;
            this.afcs = new HashMap<>(maxPeriod * maxRepeats);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(path, ex);
        }
    }

    /**
     * Dump the parameters in their text file form given the destination path.
     * @param path where to print the parameters file.
     */
    public void print(final String path) {
        try (final BufferedWriter writer = openBufferedWriter(path);
             final PrintWriter printWriter = new PrintWriter(writer)) {
            print(printWriter, path);
            this.path = path;
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(path, ex);
        }
    }

    /**
     * Dump the parameters in their text file form given the sink print writer.
     * @param printWriter the sink for the dump.
     * @param path if not null indicates the input writer resource.
     */
    public void print(final PrintWriter printWriter, final String path) {
        final StringBuilder lineBuilder = new StringBuilder(LINE_BUILDER_BUFFER_SIZE);
        lineBuilder.append(String.format("%5s", "1"));
        for (int i = 2; i <= maxRepeats; i++) {
            lineBuilder.append("  ");
            lineBuilder.append(String.format("%5s", i));
        }
        printWriter.println(lineBuilder.toString());
        printTable(printWriter, lineBuilder, GOP_TABLE_NAME, gop);
        printTable(printWriter, lineBuilder, GCP_TABLE_NAME, gcp);
        printTable(printWriter, lineBuilder, API_TABLE_NAME, api);
        if (path != null) {
            this.path = path;
        }
    }

    private void printTable(final PrintWriter printWriter, final StringBuilder lineBuilder, final String tableName, final double[][] table)  {
        printWriter.println(tableName + ":");
        for (final double[] row : table) {
            lineBuilder.setLength(0); // i.e. clear the builder.
            lineBuilder.append(String.format("%5s", String.format("%.2f", row[0])));
            for (int i = 1; i < row.length; i++) {
                lineBuilder.append("  ");
                lineBuilder.append(String.format("%5s", String.format("%.2f", row[i])));
            }
            printWriter.println(lineBuilder.toString());
        }
    }

    private DragstrParams(final int maxPeriod, final int maxRepeats, final double[][] gop, final double[][] gcp, final double[][] api) {
        this.maxPeriod = maxPeriod;
        this.maxRepeats = maxRepeats;
        this.gop = gop;
        this.gcp = gcp;
        this.api = api;
        this.afcs = new HashMap<>(maxPeriod * maxRepeats);
    }


    private static void checkMatricesAreValid(final int maxPeriod, final int maxRepeats, final double[][] gopMatrix,
                                              final double[][] gcpMatrix, final double[][] apiMatrix) {
        checkMatrixIsValid(maxPeriod, maxRepeats, gopMatrix, GOP_TABLE_NAME);
        checkMatrixIsValid(maxPeriod, maxRepeats, gcpMatrix, GCP_TABLE_NAME);
        checkMatrixIsValid(maxPeriod, maxRepeats, apiMatrix, API_TABLE_NAME);
    }

    private static void checkMatrixIsValid(final int maxPeriod, final int maxRepeatLength, final double[][] matrix, final String name) {
        Utils.nonNull(matrix, "the " + name + " matrix provided cannot be null");
        if (matrix.length != maxPeriod) {
            throw new UserException.BadInput("the " + name + " matrix provided has the wrong number of rows");
        } else {
            for (final double[] row : matrix) {
                Utils.nonNull(row, "the " + name + " matrix contains null rows");
                if (row.length != maxRepeatLength) {
                    throw new UserException.BadInput("the " + name + " matrix contains rows with length that does not match the max repeat length");
                }
            }
        }
    }

    private static double[][] mandatoryMatrix(final Map<String, double[][]> tableData, final String name) {
        final double[][] result = tableData.get(name);
        if (result == null) {
            throw new UserException.BadInput("missing matrix " + name);
        } else {
            return result;
        }
    }

    private static double[][] linesToMatrix(final List<String> lines, final int expectedNumberOfColumns) {

        final double[][] result = new double[lines.size()][expectedNumberOfColumns];
        for (int i = 0; i < lines.size(); i++) {
            final String line = lines.get(i);
            final String[] parts = line.split("\\s+");

            if (parts.length < expectedNumberOfColumns) {
                throw new UserException.BadInput("line has the wrong number of columns");
            }
            int k = 0;
            for (int j = 0; j < parts.length; j++) {
                if (!parts[j].isEmpty()) {
                    if (k >= expectedNumberOfColumns) {
                        throw new UserException.BadInput("line has the wrong number of columns");
                    }
                    try {
                        final double val = Double.parseDouble(parts[j]);
                        if (Double.isNaN(val) || Double.isInfinite(val) || val < 0.0) {
                            throw new NullPointerException();
                        }
                        result[i][k++] = val;
                    } catch (final NumberFormatException ex) {
                        throw new UserException.BadInput(
                                String.format("score is not a valid Phred value (%d,%d) == %s", i + 1, j + 1, parts[j]));
                    }
                }
            }
        }
        return result;
    }

    private double lookup(final double[][] matrix, final int period, final int repeats) {
        ParamUtils.isPositive(period, "period");
        ParamUtils.isPositive(repeats, "repeat length in units");
        final int periodIndex = period < maxPeriod ? period - 1 : maxPeriod - 1;
        final int repeatIndex = repeats < maxRepeats ? repeats - 1 : maxRepeats - 1;
        return matrix[periodIndex][repeatIndex];
    }

    /**
     * Returns the GOP for a specific period and repeat length in Phred scale.
     * @param period the query period/unit length.
     * @param repeats the query repeat length
     * @return 0 or greater.
     */
    public double gop(final int period, final int repeats) {
        return lookup(gop, period, repeats);
    }

    /**
     * Returns the GCP for a specific period and repeat length in Phred scale.
     * @param period the query period/unit length.
     * @param repeats the query repeat length
     * @return 0 or greater.
     */
    public double gcp(final int period, final int repeats) {
        return lookup(gcp, period, repeats);
    }

    /**
     * Returns the API for a specific period and repeat length in Phred scale.
     * @param period the query period/unit length.
     * @param repeats the query repeat length
     * @return 0 or greater.
     */
    public double api(final int period, final int repeats) {
        return lookup(api, period, repeats);
    }

    /**
     * Return the maximum period this parameter collection has specified values for.
     * @return 1 or greater.
     */
    public int maximumPeriod() {
        return maxPeriod;
    }

    /**
     * Return the maximum repeat length this parameter collection has specified values for.
     * @return 1 or greater.
     */
    public int maximumRepeats() {
        return maxRepeats;
    }

    public int maximumLengthInBasePairs() {
        return maxPeriod * maxRepeats;
    }

    /**
     * Return an AFC based on this DRAGstr parameters and context information on the targeted STR.
     *
     * @param period the target STR period length.
     * @param repeat the target STR repeat length.
     * @param ploidy the ploidy at that STR site.
     * @param snpHet the snp Heterozygosity at that site.
     * @param scale the scale, the weight of the prior.
     * @return never {@code null}.
     */
    public AlleleFrequencyCalculator getAFCalculator(final int period, final int repeat, final int ploidy, final double snpHet, final double scale) {
        final String keyString = "" + period + '/' + repeat + '/' + ploidy + '/' + snpHet + '/' + scale;
        return afcs.computeIfAbsent(keyString, k -> AlleleFrequencyCalculator.makeCalculator(this, period, repeat, ploidy, snpHet, scale));
    }

    @Override
    public String toString() {
        return path;
    }
}
