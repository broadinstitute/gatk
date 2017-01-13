package org.broadinstitute.hellbender.tools.coveragemodel.nd4jutils;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.IntStream;

/**
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class Nd4jIOUtils {

    /**
     * Write NDArray to a binary dump
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
     * Read NDArray from a binary dump
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
     * Write NDArray to a tab-separated file
     *
     * @param arr an instance of {@link INDArray}
     * @param outputFile the output file
     * @param columnNames column names
     * @param rowNames row names
     */
    public static void writeNDArrayToTextFile(@Nonnull final INDArray arr,
                                              @Nonnull final File outputFile,
                                              @Nullable final List<String> rowNames,
                                              @Nullable final List<String> columnNames) {
        Utils.validateArg(arr.rank() <= 2, "only rank-1 and rank-2 NDArray objects can be saved");

        final int rowDimension = arr.shape()[0];
        final int colDimension = arr.shape()[1];

        final TableColumnCollection columnNameCollection;
        final List<String> columnNameList = new ArrayList<>();
        if (columnNames == null) {
            columnNameList.add("ROW_NAME");
            IntStream.range(0, colDimension)
                    .mapToObj(colIndex -> String.format("COL_%d", colIndex))
                    .forEach(columnNameList::add);
            columnNameCollection = new TableColumnCollection(columnNameList);
        } else {
            Utils.validateArg(columnNames.size() == colDimension, "The length of column name list does not match" +
                    " the column dimension of the provided NDArray");
            columnNameList.add("ROW_NAME");
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
            for (int ri = 0; ri < rowDimension; ri++) {
                arrayWriter.writeRecord(new DoubleVectorRow(rowNamesList.get(ri), arr.getRow(ri).dup().data().asDouble()));
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "Could not write the array");
        }
    }

    /**
     * Read NDArray from a tab-separated file
     *
     * @param inputFile the input tab-separated file
     * @return an {@link INDArray}
     */
    public static INDArray readNDArrayFromTextFile(@Nonnull final File inputFile) {
        final List<INDArray> rows = new LinkedList<>();
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
            arrayReader.iterator().forEachRemaining(rowData -> rows.add(Nd4j.create(rowData.data,
                    new int[] {1, colDimension})));
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile(inputFile, "Could not read NDArray tsv file");
        }
        return Nd4j.vstack(rows);
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

        public void composeDataLine(final DataLine dataLine) {
            dataLine.append(rowName);
            dataLine.append(data);
        }
    }

}
