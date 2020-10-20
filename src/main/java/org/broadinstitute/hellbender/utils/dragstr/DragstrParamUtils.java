package org.broadinstitute.hellbender.utils.dragstr;

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.functional.IntToDoubleBiFunction;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Utils to read and write {@link DragstrParams} instances from an to files and other resources.
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
 *     Currently there is no fix/standard set of annotations but is just a mechanism to add some additional
 *     metadata that might be useful for tracking purposes.
 * </p>
 */
public final class DragstrParamUtils {

    public static final String GOP_TABLE_NAME = "GOP";
    public static final String GCP_TABLE_NAME = "GCP";
    public static final String API_TABLE_NAME = "API";

    /**
     * Initial size of the buffer/builder used to compose lines in the output text format.
     */
    private static int LINE_BUILDER_BUFFER_SIZE = 1024;


    public static DragstrParams parse(final GATKPath path) {
        if (path == null) {
            return null;
        }
        try (final BufferedReader reader = openBufferedReader(path.toString())) {
            return parse(reader, path.toString());
        } catch (final IOException e) {
            throw  new UserException.CouldNotReadInputFile(path.toString(), e);
        }
    }

    public static void print(final DragstrParams params, final GATKPath path, final Object ... annotations) {
        Utils.nonNull(params, "the input params cannot be null");
        Utils.nonNull(path, "the input path cannot be null");
        try (final PrintWriter pw = new PrintWriter(path.getOutputStream())) {
            print(params, pw, annotations);
        }
        params.setName(path.toString());
    }

    private static void print(final DragstrParams params, final PrintWriter writer, final Object ... annotations) {
        writer.println("############################################################################################");
        writer.println("# DragstrParams");
        writer.println("# -------------------------");
        for (int i = 0; i < annotations.length;) {
            final Object name = annotations[i++];
            final Object value = i < annotations.length ? annotations[i++] : null;
            writer.println("# " + name + (value != null ? (" = " + value) : ""));
        }
        writer.println("############################################################################################");
        printTables(params, writer);
    }

    private static BufferedReader openBufferedReader(final String path) {
        try {
            return Files.newBufferedReader(Paths.get(path));
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(path, ex);
        }
    }

    private static DragstrParams parse(final BufferedReader reader, final String name) {
        try {
            String header;
            while ((header = reader.readLine()) != null) {
                if (!header.startsWith("#")) {
                    break;
                }
            }
            if (header == null) {
                throw new UserException.BadInput("there is no content in the dragstr-params file " + name);
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
            return DragstrParams.of(maxPeriod, maxRepeats, gopMatrix, gcpMatrix, apiMatrix, name);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(name, ex);
        }
    }

    private static void printTable(final PrintWriter printWriter, final StringBuilder lineBuilder, final String tableName,
                            final int maxPeriod, final int maxRepeat, IntToDoubleBiFunction func)  {
        printWriter.println(tableName + ":");
        for (int i = 1; i <= maxPeriod; i++) {
            lineBuilder.setLength(0);
            lineBuilder.append(String.format("%5s", String.format("%.2f", func.apply(i, 1))));
            for (int j = 2; j <= maxRepeat; j++) {
                lineBuilder.append("  ");
                lineBuilder.append(String.format("%5s", String.format("%.2f", func.apply(i, j))));
            }
            printWriter.println(lineBuilder.toString());
        }
    }

    /**
     * Dump the parameters in their text file form given the sink print writer.
     * @param printWriter the sink for the dump.
     *
     */
    private static void printTables(final DragstrParams params, final PrintWriter printWriter) {
        final StringBuilder lineBuilder = new StringBuilder(LINE_BUILDER_BUFFER_SIZE);
        lineBuilder.append(String.format("%5s", "1"));
        for (int i = 2; i <= params.maximumRepeats(); i++) {
            lineBuilder.append("  ");
            lineBuilder.append(String.format("%5s", i));
        }
        printWriter.println(lineBuilder.toString());
        printTable(printWriter, lineBuilder, GOP_TABLE_NAME, params.maximumPeriod(), params.maximumRepeats(), params::gop);
        printTable(printWriter, lineBuilder, GCP_TABLE_NAME, params.maximumPeriod(), params.maximumRepeats(), params::gcp);
        printTable(printWriter, lineBuilder, API_TABLE_NAME, params.maximumPeriod(), params.maximumRepeats(), params::api);
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
}
