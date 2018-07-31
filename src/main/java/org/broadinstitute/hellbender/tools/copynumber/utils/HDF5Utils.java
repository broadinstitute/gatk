package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * TODO move into hdf5-java-bindings
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HDF5Utils {
    private static final Logger logger = LogManager.getLogger(HDF5Utils.class);

    //writing intervals as a string matrix is expensive,
    //so we instead store a map from integer indices to contig strings and
    //store (index, start, end) in a double matrix;
    //these are stored in the following sub-paths
    public static final String INTERVAL_CONTIG_NAMES_SUB_PATH = "/indexed_contig_names";
    public static final String INTERVAL_MATRIX_SUB_PATH = "/transposed_index_start_end";

    //Java HDF5 has a hard limit on the number of elements in a matrix given by the following
    public static final int MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX = Integer.MAX_VALUE / Byte.SIZE;
    //this limit requires us to write matrices exceeding the limit as sets of submatrix chunks;
    //we use the following sub-paths to store the necessary information
    public static final String NUMBER_OF_ROWS_SUB_PATH = "/num_rows";
    public static final String NUMBER_OF_COLUMNS_SUB_PATH = "/num_columns";
    public static final String NUMBER_OF_CHUNKS_SUB_PATH = "/num_chunks";
    public static final String CHUNK_INDEX_PATH_SUFFIX = "/chunk_";

    private enum IntervalField {
        CONTIG_INDEX (0),
        START (1),
        END (2);

        private final int index;

        IntervalField(final int index) {
            this.index = index;
        }
    }
    private static final int NUM_INTERVAL_FIELDS = IntervalField.values().length;

    private HDF5Utils() {}

    /**
     * Reads a list of intervals from an HDF5 file using the sub-paths and conventions used by {@link #writeIntervals}.
     */
    public static List<SimpleInterval> readIntervals(final HDF5File file,
                                                     final String path) {
        final String[] contigNames = file.readStringArray(path + INTERVAL_CONTIG_NAMES_SUB_PATH);
        final double[][] matrix = file.readDoubleMatrix(path + INTERVAL_MATRIX_SUB_PATH);
        final int numIntervals = matrix[0].length;
        return IntStream.range(0, numIntervals)
                .mapToObj(i -> (new SimpleInterval(
                        contigNames[(int) matrix[IntervalField.CONTIG_INDEX.index][i]],
                        (int) matrix[IntervalField.START.index][i],
                        (int) matrix[IntervalField.END.index][i])))
                .collect(Collectors.toList());
    }

    /**
     * Given an HDF5 file and an HDF5 path, writes a list of intervals to hard-coded sub-paths.
     * Contig names are represented by a string array, while intervals are represented by a double matrix,
     * in which the contigs are represented by their index in the aforementioned string array.
     */
    public static <T extends SimpleInterval> void writeIntervals(final HDF5File file,
                                                                 final String path,
                                                                 final List<T> intervals) {
        final Map<String, Double> contigNamesToIndexMap = new LinkedHashMap<>();
        final double[][] matrix = new double[NUM_INTERVAL_FIELDS][intervals.size()];
        for (int i = 0; i < intervals.size(); i++) {
            final SimpleInterval interval = intervals.get(i);
            contigNamesToIndexMap.putIfAbsent(interval.getContig(), (double) contigNamesToIndexMap.keySet().size());
            matrix[IntervalField.CONTIG_INDEX.index][i] = contigNamesToIndexMap.get(interval.getContig());
            matrix[IntervalField.START.index][i] = interval.getStart();
            matrix[IntervalField.END.index][i] = interval.getEnd();
        }
        file.makeDoubleMatrix(path + INTERVAL_MATRIX_SUB_PATH, matrix);
        file.makeStringArray(path + INTERVAL_CONTIG_NAMES_SUB_PATH,
                contigNamesToIndexMap.keySet().toArray(new String[contigNamesToIndexMap.keySet().size()]));
    }

    /**
     * Reads a large matrix stored as a set of chunks (submatrices) using the sub-paths and conventions
     * used by {@link #writeChunkedDoubleMatrix}.
     */
    public static double[][] readChunkedDoubleMatrix(final HDF5File file,
                                                     final String path) {
        Utils.nonNull(file);
        IOUtils.canReadFile(file.getFile());
        Utils.nonNull(path);

        final String numRowsPath = path + NUMBER_OF_ROWS_SUB_PATH;
        final String numColumnsPath = path + NUMBER_OF_COLUMNS_SUB_PATH;
        final String numChunksPath = path + NUMBER_OF_CHUNKS_SUB_PATH;
        Utils.validateArg(file.isPresent(numRowsPath) && file.isPresent(numColumnsPath) && file.isPresent(numChunksPath),
                String.format("HDF5 file %s does not contain a chunked matrix in path %s.", file.getFile().getAbsolutePath(), path));

        final int numRows = (int) file.readDouble(numRowsPath);
        final int numColumns = (int) file.readDouble(numColumnsPath);
        final int numChunks = (int) file.readDouble(numChunksPath);

        final double[][] fullMatrix = new double[numRows][numColumns];
        int numRowsRead = 0;
        for (int chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
            final double[][] matrixChunk = file.readDoubleMatrix(path + CHUNK_INDEX_PATH_SUFFIX + chunkIndex);
            if (numRowsRead + matrixChunk.length > numRows) {
                throw new UserException.BadInput("Matrix chunk contains too many rows.");
            }
            if (matrixChunk[0].length != numColumns) {
                throw new UserException.BadInput("Matrix chunk does not contain expected number of columns.");
            }
            System.arraycopy(matrixChunk, 0, fullMatrix, numRowsRead, matrixChunk.length);
            numRowsRead += matrixChunk.length;
        }
        if (numRowsRead != numRows) {
            throw new UserException.BadInput("Matrix chunks do not contain expected total number of rows.");
        }
        return fullMatrix;
    }

    /**
     * Given a large matrix, chunks the matrix into equally sized subsets of rows
     * (plus a subset containing the remainder, if necessary) and writes these submatrices to indexed sub-paths
     * to avoid a hard limit in Java HDF5 on the number of elements in a matrix given by
     * {@code MAX_NUM_VALUES_PER_HDF5_MATRIX}. The number of chunks is determined by {@code maxChunkSize},
     * which should be set appropriately for the desired number of columns.
     *
     * @param maxChunkSize  The maximum number of values in each chunk. Decreasing this number will reduce
     *                      heap usage when writing chunks, which requires subarrays to be copied.
     *                      However, since a single row is not allowed to be split across multiple chunks,
     *                      the number of columns must be less than the maximum number of values in each chunk.
     */
    public static void writeChunkedDoubleMatrix(final HDF5File file,
                                                final String path,
                                                final double[][] matrix,
                                                final int maxChunkSize) {
        Utils.nonNull(file);
        IOUtils.canReadFile(file.getFile());
        Utils.nonNull(path);
        Utils.nonNull(matrix);
        ParamUtils.inRange(maxChunkSize, 1 , MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX,
                String.format("Maximum chunk size must be in [1, %d].", MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX));
        final long numRows = matrix.length;
        Utils.validateArg(numRows > 0, "Matrix must contain at least one row.");
        final long numColumns = matrix[0].length;
        Utils.validateArg(numColumns > 0, "Matrix must contain at least one column.");
        Utils.validateArg(numColumns <= maxChunkSize,
                String.format("Number of columns (%d) exceeds the maximum number of values allowed per chunk (%d).",
                        numColumns, maxChunkSize));

        final int numRowsPerFilledChunk = (int) (maxChunkSize / numColumns);
        final int numFilledChunks = numRowsPerFilledChunk == 0 ? 0 : (int) numRows / numRowsPerFilledChunk;
        final boolean needPartialChunk = numFilledChunks == 0 || numRows % numRowsPerFilledChunk != 0;

        logger.debug("Number of values in matrix / maximum number allowed for HDF5 matrix: " + (double) numRows * numColumns / MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX);
        logger.debug("Maximum number of values per chunk: " + maxChunkSize);
        logger.debug("Number of filled chunks: " + numFilledChunks);
        logger.debug("Number of rows per filled chunk: " + numRowsPerFilledChunk);
        logger.debug("Partial chunk needed: " + needPartialChunk);

        final String numRowsPath = path + NUMBER_OF_ROWS_SUB_PATH;
        final String numColumnsPath = path + NUMBER_OF_COLUMNS_SUB_PATH;
        final String numChunksPath = path + NUMBER_OF_CHUNKS_SUB_PATH;
        file.makeDouble(numRowsPath, numRows);
        file.makeDouble(numColumnsPath, numColumns);
        file.makeDouble(numChunksPath, needPartialChunk ? numFilledChunks + 1 : numFilledChunks);

        //TODO we could add makeDoubleMatrix(path, matrix, startRow, endRow, startCol, endCol) method to avoid copying
        int numRowsWritten = 0;
        for (int chunkIndex = 0; chunkIndex < numFilledChunks; chunkIndex++) {
            final double[][] matrixChunk = new double[numRowsPerFilledChunk][(int) numColumns];
            System.arraycopy(matrix, numRowsWritten, matrixChunk, 0, numRowsPerFilledChunk);
            file.makeDoubleMatrix(path + CHUNK_INDEX_PATH_SUFFIX + chunkIndex, matrixChunk);    //write filled chunks
            numRowsWritten += numRowsPerFilledChunk;
        }
        if (needPartialChunk) {
            final int numRowsPartialChunk = (int) numRows - numRowsWritten;
            logger.debug("Number of rows in partial chunk: " + numRowsPartialChunk);
            final double[][] matrixChunk = new double[numRowsPartialChunk][(int) numColumns];
            System.arraycopy(matrix, numRowsWritten, matrixChunk, 0, numRowsPartialChunk);
            file.makeDoubleMatrix(path + CHUNK_INDEX_PATH_SUFFIX + numFilledChunks, matrixChunk);    //write final partially filled chunk
        }
    }
}
