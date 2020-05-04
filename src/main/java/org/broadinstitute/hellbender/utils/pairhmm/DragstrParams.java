package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.exceptions.UserException;
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

public final class DragstrParams {

    public static final DragstrParams DEFAULT = new DragstrParams(
            DragstrConstants.DEFAULT_MAX_PERIOD,
            DragstrConstants.DEFAULT_MAX_REPEATS,
            new double[][] { /* GOP */
                    {45.00, 45.00, 45.00, 45.00, 45.00, 45.00, 40.50, 33.50, 28.00, 24.00, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75},
                    {39.50, 39.50, 39.50, 39.50, 36.00, 30.00, 27.25, 25.00, 24.25, 24.75, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.75},
                    {38.50, 41.00, 41.00, 41.00, 41.00, 37.50, 35.25, 34.75, 34.75, 33.25, 33.25, 33.25, 32.50, 30.75, 28.50, 29.00, 29.00, 29.00, 29.00, 29.00},
                    {37.50, 39.00, 39.00, 37.75, 34.00, 34.00, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 31.75, 31.75, 31.75, 31.75, 31.75},
                    {37.00, 40.00, 40.00, 40.00, 36.00, 35.00, 24.50, 24.50, 24.50, 24.50, 22.50, 22.50, 22.50, 23.50, 23.50, 23.50, 23.50, 23.50, 23.50, 23.50},
                    {36.25, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00},
                    {36.00, 40.50, 40.50, 40.50, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75},
                    {36.25, 39.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75}},
            new double[][] { /* GCP */
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 10.0),
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 5.0),        // = 10.0 / 2
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 10.0 / 3.0),
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 10.0 / 4.0),
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 2.0),        // = 10.0 / 5
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 10.0 / 6.0),
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 10.0 / 7.0),
                    MathUtils.doubles(DragstrConstants.DEFAULT_MAX_REPEATS, 10.0 / 8.0),
            },
            new double[][] { /* API */
                    {39.00, 39.00, 37.00, 35.00, 32.00, 26.00, 20.00, 16.00, 12.00, 10.00, 8.00, 7.00, 7.00, 6.00, 6.00, 5.00, 5.00, 4.00, 4.00, 4.00},
                    {30.00, 30.00, 29.00, 22.00, 17.00, 14.00, 11.00, 8.00, 6.00, 5.00, 4.00, 4.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 2.00, 2.00},
                    {27.00, 27.00, 25.00, 18.00, 14.00, 12.00, 9.00, 7.00, 5.00, 4.00, 3.00, 3.00, 3.00, 3.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00},
                    {27.00, 27.00, 18.00, 9.00, 9.00, 9.00, 9.00, 3.00, 3.00, 3.00, 3.00, 3.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00},
                    {29.00, 29.00, 18.00, 8.00, 8.00, 8.00, 4.00, 3.00, 3.00, 3.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00},
                    {25.00, 25.00, 10.00, 10.00, 10.00, 4.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00},
                    {21.00, 21.00, 11.00, 11.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00},
                    {18.00, 18.00, 10.00, 6.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00}}
    );


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

    public static DragstrParams of(int maxPeriod, int maxRepeats, double[][] gop, double[][] gcp, double[][] api) {
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
            final String header = reader.readLine();
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
            if (tableName == null) {
                throw new UserException.BadInput("table with no name");
            }
            tables.put(tableName, linesToMatrix(tableLines, repeats.length));
            final double[][] gopMatrix = mandatoryMatrix(tables, "GOP");
            final double[][] gcpMatrix = mandatoryMatrix(tables, "GCP");
            final double[][] apiMatrix = mandatoryMatrix(tables, "API");
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

    public void print(final String path) {
        try (final BufferedWriter writer = openBufferedWriter(path);
             final PrintWriter printWriter = new PrintWriter(writer)) {
            final StringBuilder lineBuilder = new StringBuilder(1024);
            lineBuilder.append(String.format("%5s", "1"));
            for (int i = 2; i <= maxRepeats; i++) {
                lineBuilder.append("  ");
                lineBuilder.append(String.format("%5s", i));
            }
            printWriter.println(lineBuilder.toString());
            printTable(printWriter, lineBuilder, "GOP", gop);
            printTable(printWriter, lineBuilder, "GCP", gcp);
            printTable(printWriter, lineBuilder, "API", api);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(path, ex);
        }
    }

    private void printTable(final PrintWriter printWriter, final StringBuilder lineBuilder, final String tableName, final double[][] table)  {
        printWriter.println(tableName + ":");
        for (final double[] row : table) {
            lineBuilder.setLength(0);
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
        if (gopMatrix.length != maxPeriod) {
            throw new IllegalArgumentException();
        } else if (gcpMatrix.length != maxPeriod) {
            throw new UserException.BadInput("the GCP matrix must have the same number of rows as the GOP");
        } else if (apiMatrix.length != maxPeriod) {
            throw new UserException.BadInput("the API matrix must have the same number of rows as the GOP");
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

    public double gop(final int period, final int repeats) {
        return lookup(gop, period, repeats);
    }

    public double gcp(final int period, final int repeats) {
        return lookup(gcp, period, repeats);
    }

    public double api(final int period, final int repeats) {
        return lookup(api, period, repeats);
    }

    public int maximumPeriod() {
        return maxPeriod;
    }

    public int maximumRepeats() {
        return maxRepeats;
    }

    public AlleleFrequencyCalculator getAFCalculator(final int period, final int repeat, final int ploidy, final double snpHet, final double scale) {
        final String keyString = "" + period + '/' + repeat + '/' + ploidy + '/' + snpHet;
        return afcs.computeIfAbsent(keyString, k -> AlleleFrequencyCalculator.makeCalculator(this, period, repeat, ploidy, snpHet, scale));
    }
}
