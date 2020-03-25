package org.broadinstitute.hellbender.utils.dragstr;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.dragstr.DragstrHyperParameters;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.stream.IntStream;

/**
 * Holds the values of the DRAGstr model parameters for different combinations of repeat unit length (period) and
 * number of repeat units.
 */
public final class DragstrParams {


    /**
     * String returned by {@link #toString()} when this object has not been retrieved from nor persisted into a file.
     */
    private static final String NON_PERSISTENT_NAME = "<non-persistent>";

    /**
     * String returned by {@link #toString()} for the default model parameters
     */
    private static final String DEFAULT_NAME = "<default>";

     /**
     * Holds the path of the file this param where loaded from or persistent into last.
     */
    private String name;

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
                   DEFAULT_API,
            DEFAULT_NAME);

    private final int maxPeriod;
    private final int maxRepeats;
    private final double[][] gop;
    private final double[][] gcp;
    private final double[][] api;

    public static DragstrParams of(final int maxPeriod, final int maxRepeats, final double[][] gop, final double[][] gcp, final double[][] api) {
        return of(maxPeriod, maxRepeats, gop, gcp, api, NON_PERSISTENT_NAME);
    }

    public static DragstrParams of(final int maxPeriod, final int maxRepeats, final double[][] gop, final double[][] gcp, final double[][] api, final String name) {
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
        checkMatricesAreValid(maxPeriod, maxRepeats, gop, gcp, api);
        return new DragstrParams(maxPeriod, maxRepeats, gop.clone(), gcp.clone(), api.clone(), name);
    }

    private DragstrParams(final int maxPeriod, final int maxRepeats, final double[][] gop, final double[][] gcp, final double[][] api, final String name) {
        this.maxPeriod = maxPeriod;
        this.maxRepeats = maxRepeats;
        this.gop = gop;
        this.gcp = gcp;
        this.api = api;
        this.name = name;
    }


    private static void checkMatricesAreValid(final int maxPeriod, final int maxRepeats, final double[][] gopMatrix,
                                              final double[][] gcpMatrix, final double[][] apiMatrix) {
        checkMatrixIsValid(maxPeriod, maxRepeats, gopMatrix, DragstrParamUtils.GOP_TABLE_NAME);
        checkMatrixIsValid(maxPeriod, maxRepeats, gcpMatrix, DragstrParamUtils.GCP_TABLE_NAME);
        checkMatrixIsValid(maxPeriod, maxRepeats, apiMatrix, DragstrParamUtils.API_TABLE_NAME);
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

    @Override
    public String toString() {
        return name;
    }

    void setName(final String name) {
        this.name = Utils.nonNull(name);
    }
}
