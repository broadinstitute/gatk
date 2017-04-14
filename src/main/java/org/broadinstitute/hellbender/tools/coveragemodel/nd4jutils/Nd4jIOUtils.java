package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import avro.shaded.com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.ImmutableTriple;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.*;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class Nd4jIOUtils {

    /**
     * Writes NDArray to a binary dump.
     *
     * @param arr an instance of {@link INDArray}
     * @param outputFile the output path
     */
    public static void writeNDArrayToBinaryDumpFile(@Nonnull final INDArray arr,
                                                    @Nonnull final File outputFile) {
        try (final DataOutputStream dos = new DataOutputStream(new FileOutputStream(outputFile))) {
            Nd4j.write(arr, dos);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "Could not write the NDArray");
        }
    }

    /**
     * Reads NDArray from a binary dump.
     *
     * @param inputFile the input Nd4j binary dump file
     * @return an {@link INDArray}
     */
    public static INDArray readNDArrayFromBinaryDumpFile(@Nonnull final File inputFile) {
        try (final DataInputStream dis = new DataInputStream(new FileInputStream(inputFile))) {
            return Nd4j.read(dis);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, "Could not read the input file");
        }
    }

    /**
     * Writes NDArray to a tab-separated file.
     *
     * @param arr an instance of {@link INDArray}
     * @param outputFile the output file
     * @param identifier the identifier string (just a name)
     * @param columnNames column names
     * @param rowNames row names
     */
    public static void writeNDArrayMatrixToTextFile(@Nonnull final INDArray arr,
                                                    @Nonnull final File outputFile,
                                                    @Nonnull final String identifier,
                                                    @Nullable final List<String> rowNames,
                                                    @Nullable final List<String> columnNames,
                                                    @Nullable final int[] originalShape) {
        Utils.nonNull(arr, "The NDArray to be written to file must be non-null");
        Utils.nonNull(outputFile, "The output file must be non-null");
        Utils.nonNull(identifier, "The identifier must be non-null");
        Utils.validateArg(arr.rank() <= 2, "only rank-1 and rank-2 NDArray objects can be saved");
        final int[] shape;
        if (originalShape != null) {
            Utils.validateArg(Arrays.stream(originalShape).reduce(1, (a, b) -> a*b) ==
                    Arrays.stream(arr.shape()).reduce(1, (a, b) -> a*b), "The given shape does not match the number" +
                    " of elements in the given array");
            shape = originalShape;
        } else {
            shape = arr.shape();
        }

        final int rowDimension = arr.shape()[0];
        final int colDimension = arr.shape()[1];

        final TableColumnCollection columnNameCollection;
        final List<String> columnNameList = new ArrayList<>();
        if (columnNames == null) {
            columnNameList.add(identifier);
            IntStream.range(0, colDimension)
                    .mapToObj(colIndex -> String.format("COL_%d", colIndex))
                    .forEach(columnNameList::add);
            columnNameCollection = new TableColumnCollection(columnNameList);
        } else {
            Utils.validateArg(columnNames.size() == colDimension, "The length of column name list does not match" +
                    " the column dimension of the provided NDArray");
            columnNameList.add(identifier);
            columnNameList.addAll(columnNames);
            columnNameCollection = new TableColumnCollection(columnNameList);
        }

        final List<String> rowNamesList = new ArrayList<>();
        if (rowNames == null) {
            IntStream.range(0, rowDimension)
                    .mapToObj(colIndex -> String.format("ROW_%d", colIndex))
                    .forEach(rowNamesList::add);
        } else {
            Utils.validateArg(rowNames.size() == rowDimension, "The length of row name string array does not match" +
                    " the column dimension of the provided NDArray");
            rowNamesList.addAll(rowNames);
        }

        try (final TableWriter<DoubleVectorRow> arrayWriter = new TableWriter<DoubleVectorRow>(outputFile, columnNameCollection) {
            @Override
            protected void composeLine(DoubleVectorRow record, DataLine dataLine) {
                record.composeDataLine(dataLine);
            }}) {
            /* write the shape info as a comment */
            arrayWriter.writeComment(getStringReprFromIntArray(shape));
            /* write rows */
            for (int ri = 0; ri < rowDimension; ri++) {
                arrayWriter.writeRecord(new DoubleVectorRow(rowNamesList.get(ri), arr.getRow(ri).dup().data().asDouble()));
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "Could not write the array");
        }
    }

    public static void writeNDArrayMatrixToTextFile(@Nonnull final INDArray arr,
                                                    @Nonnull final File outputFile,
                                                    @Nonnull final String identifier,
                                                    @Nullable final List<String> rowNames,
                                                    @Nullable final List<String> columnNames) {
        writeNDArrayMatrixToTextFile(arr, outputFile, identifier, rowNames, columnNames, null);
    }

    /**
     * Reads NDArray from a tab-separated file and returns a triple of the array, row names, and column names
     *
     * @param inputFile the input tab-separated file
     * @return a triple of ({@link INDArray}, row names, column names)
     */
    public static ImmutableTriple<INDArray, List<String>, List<String>> readNDArrayMatrixFromTextFileWithRowAndColumnNames(
            @Nonnull final File inputFile) {
        final List<INDArray> rows = new LinkedList<>();
        final List<String> rowNames = new ArrayList<>();
        final List<String> columnNames;
        try (final TableReader<DoubleVectorRow> arrayReader = new TableReader<DoubleVectorRow>(inputFile) {
            @Override
            protected DoubleVectorRow createRecord(DataLine dataLine) {
                final int colDimension = columns().columnCount() - 1;
                final String[] dataLineToString = dataLine.toArray();
                if (dataLineToString.length != colDimension + 1) {
                    throw new UserException.BadInput("The input NDArray tsv file is malformed");
                } else {
                    final double[] rowData = Arrays.stream(dataLineToString, 1, colDimension + 1)
                            .mapToDouble(Double::new).toArray();
                    return new DoubleVectorRow(dataLineToString[0], rowData);
                }
            }
        }) {
            final int colDimension = arrayReader.columns().columnCount() - 1;
            columnNames = arrayReader.columns().names();
            arrayReader.iterator().forEachRemaining(rowData -> {
                rows.add(Nd4j.create(rowData.data, new int[]{1, colDimension}));
                rowNames.add(rowData.rowName);
            });
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, "Could not read NDArray tsv file");
        }
        return ImmutableTriple.of(Nd4j.vstack(rows), rowNames, columnNames);
    }

    /**
     * Reads NDArray from a tab-separated file and returns an instance of {@link INDArray}.
     *
     * @param inputFile the input tab-separated file
     * @return an {@link INDArray}
     */
    public static INDArray readNDArrayMatrixFromTextFile(@Nonnull final File inputFile) {
        return readNDArrayMatrixFromTextFileWithRowAndColumnNames(inputFile).getLeft();
    }

    /**
     * Writes a tensor NDArray (rank >= 2) to a tab-separated file.
     *
     * For a tensor of dim D, it is flattened along the last D - 1 dimensions (in c order) and is written as a
     * matrix. The shape of the tensor is written as a comment line for proper reshaping upon loading.
     *
     * @param arr an arbitrary NDArray
     * @param outputFile output file
     * @param identifier an identifier string
     * @param rowNames list of row names
     */
    public static void writeNDArrayTensorToTextFile(@Nonnull final INDArray arr,
                                                    @Nonnull final File outputFile,
                                                    @Nonnull final String identifier,
                                                    @Nullable final List<String> rowNames) {
        Utils.nonNull(arr, "The NDArray to be written to file must be non-null");
        Utils.nonNull(outputFile, "The output file must be non-null");
        Utils.validateArg(arr.rank() >= 2, "The array must be at least rank-2");

        final int[] shape = arr.shape();
        int extraDims = 1;
        for (int i = 1; i < arr.rank(); i++) {
            extraDims *= shape[i];
        }
        final INDArray tensorToMatrix = arr.reshape('c', new int[] {shape[0], extraDims});
        writeNDArrayMatrixToTextFile(tensorToMatrix, outputFile, identifier, rowNames, null, arr.shape());
    }

    /**
     * Reads an NDArray of rank >= 2 from a tab-separated file.
     *
     * @param inputFile the input tab-separated file
     * @return an instance of {@link INDArray}
     */
    public static INDArray readNDArrayTensorFromTextFile(@Nonnull final File inputFile) {
        Utils.regularReadableUserFile(inputFile);

        /* read the first line, make sure it is a header line, and get the dimension of the tensor */
        final String commentLine;
        try (final BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
            commentLine = reader.readLine();
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read the input file");
        }
        Utils.validateArg(commentLine.startsWith(TableUtils.COMMENT_PREFIX),
                "Invalid input file: the first line must be a comment line (starting with " + TableUtils.COMMENT_PREFIX
                + ") followed by the tensor shape");
        final int shape[] = getIntArrayFromStringRepr(commentLine.substring(TableUtils.COMMENT_PREFIX.length()));
        final INDArray tensorToMatrix = readNDArrayMatrixFromTextFile(inputFile);
        return tensorToMatrix.reshape('c', shape);
    }

    /**
     * Takes an int array and forms a string identifier out of it
     *
     * @param indices an array of indices
     * @return a string representation
     */
    @VisibleForTesting
    static String getStringReprFromIntArray(int... indices) {
        Utils.nonNull(indices, "The indices array must be non-null");
        Utils.validateArg(indices.length > 0, "The indices array must be non-empty");
        return Arrays.stream(indices).mapToObj(Integer::toString).collect(Collectors.joining("_", "[", "]"));
    }

    /**
     * Takes a string representation of an int array and returns an int array
     *
     * @param stringRepr the string representation
     * @return an int array
     */
    @VisibleForTesting
    static int[] getIntArrayFromStringRepr(@Nonnull final String stringRepr) {
        Utils.nonNull(stringRepr, "The column name must be non-null");
        Utils.validateArg(stringRepr.length() > 2 && stringRepr.startsWith("[") && stringRepr.endsWith("]"),
                "The provided string is not a valid shape string");
        return Arrays.stream(stringRepr.substring(1, stringRepr.length() - 1).split("_"))
                .mapToInt(Integer::valueOf).toArray();
    }

    /**
     * A named double array, represents a row in a tab-separated text representation of a matrix
     */
    private static final class DoubleVectorRow {
        final String rowName;
        final double[] data;

        DoubleVectorRow(final String rowName, final double[] data) {
            this.rowName = rowName;
            this.data = data;
        }

        void composeDataLine(final DataLine dataLine) {
            dataLine.append(rowName);
            dataLine.append(data);
        }
    }

}
