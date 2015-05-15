package org.broadinstitute.hellbender.utils.hdf5;

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.HDF5Constants;
import ncsa.hdf.hdf5lib.exceptions.HDF5LibraryException;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.File;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Stream;

/**
 * HDF5 File reader.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5Reader implements AutoCloseable {

    /**
     * Reference to the underlying HDF5 file location.
     */
    private final File file;

    /**
     * The underlying HDF5 file open resource id as provided by {@link H5#H5Fopen}.
     *
     * <p>
     *     A value of {@code -1} indicates that the file is closed.
     * </p>
     */
    private int fileId;

    /**
     * Special {@link #fileId} value to indicate that the reader is closed.
     */
    private int FILE_ID_WHEN_CLOSED = -1;

    /**
     * Creates a new HDF5 reader on a existing file.
     *
     * @param file the target file.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws GATKException if the HDF5 library is not supported or could not be initialized.
     */
    public HDF5Reader(final File file) {
        HDF5Library.getLibrary();
        if (file == null) {
            throw new IllegalArgumentException("the file cannot be null.");
        }
        this.file = file;
        try {
            fileId = H5.H5Fopen(file.getAbsolutePath(), HDF5Constants.H5F_ACC_RDONLY, HDF5Constants.H5P_DEFAULT);
        } catch (final HDF5LibraryException e) {
            throw new GATKException(
                    String.format("exception when opening '%s' for read-only access: %s",file.getAbsolutePath(), e.getMessage()),e);
        }
        if (fileId < 0) {
            throw new GATKException(
                    String.format("failure when opening '%s' for read-only access; negative fileId: %d",file.getAbsolutePath(),fileId)
            );
        }
    }

    /**
     * Returns a reference to the underlying HDF5 file location.
     * @return never {@code null}.
     */
    public File getFile() {
        return file;
    }

    /**
     * Close this file reader.
     *
     * <p>
     *     Further file read operations will result in a {@link IllegalStateException}.
     * </p>
     */
    public void close() {
        if (isClosed()) {
            return;
        }
        try {
            H5.H5Fclose(fileId);
            fileId = FILE_ID_WHEN_CLOSED;
        } catch (final HDF5LibraryException e) {
            throw new GATKException(
                    String.format("failure when closing '%s' from read-only access: %s",file.getAbsolutePath(),e.getMessage()),e);
        }
    }

    /**
     * Checks whether the reader has been closed.
     * @return {@code true} iff the read is closed.
     */
    protected boolean isClosed() {
        return fileId == FILE_ID_WHEN_CLOSED;
    }

    /**
     * Retrieves a string array given its full-path within the HDF5 file.
     * @param fullPath the target data-set name.
     * @throws GATKException if the operation failed, for example if there is no such a data-set
     * @return never {@code null}.
     */
    public String[] readStringArray(final String fullPath) {
        return readDataset(fullPath, (dataSetId, typeId, dimensions) -> {
            if (dimensions.length != 1) {
                throw new GATKException(
                        String.format("expected 1-D array for data-set '%s' in '%s' but it is %d-D", fullPath, file, dimensions.length));
            }
            final String[] result = new String[(int) dimensions[0]];
            final int code = H5.H5Dread_string(dataSetId, typeId, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, result);
            if (code < 0) {
                throw new GATKException(String.format("getting strings from data-set '%s' in file '%s' resulted in code: %d", fullPath, file, code));
            }
            return result;
        });
    }
    /**
     * Reads a double value from a particular position in the underlying HDF5 file.
     *
     * @param fullPath the path to the double value in the HDF5 file.
     * @return the stored value.
     * @throws IllegalArgumentException if {@code fullPath} is {@code null}.
     * @throws GATKException if {@code fullPath} does not exist, contains the wrong data type (non-double) or
     *    is an array with several values or multidimensional.
     */
    public double readDouble(final String fullPath) {
        return readDataset(fullPath, (dataSetId, typeId, dimensions) -> {
            if (dimensions.length != 1) {
                throw new GATKException(
                        String.format("expected 1-D array for data-set '%s' in '%s' but it is %d-D", fullPath, file, dimensions.length));
            }
            if (dimensions[0] != 1) {
                throw new GATKException(
                        String.format("expected single value array for data-set '%s' in '%s' but it has %d values", fullPath, file, dimensions[0]));
            }
            final double[] values = new double[1];
            final int code = H5.H5Dread_double(dataSetId, typeId, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, values);
            if (code < 0) {
                throw new GATKException(String.format("getting a double from data-set '%s' in file '%s' resulted in code: %d", fullPath, file, code));
            }
            return values[0];
        });
    }

    /**
     * Reads a double array from a particular position in the underlying HDF5 file.
     *
     * @param fullPath the path.
     * @return never {@code null}, non-existing data or the wrong type in the HDF5
     * file will result in a {@link GATKException} instead.
     * @throws IllegalArgumentException if {@code fullPath} is {@code null}.
     * @throws GATKException if {@code fullPath} does not exist, contains the wrong data type (non-double) or
     *    is multidimensional.
     */
    public double[] readDoubleArray(final String fullPath) {
        return readDataset(fullPath, (dataSetId, typeId, dimensions) -> {
            if (dimensions.length != 1) {
                throw new GATKException(
                        String.format("expected 1-D array for data-set '%s' in '%s' but it is %d-D", fullPath, file, dimensions.length));
            }
            final double[] result = new double[(int) dimensions[0]];
            final int code = H5.H5Dread_double(dataSetId, typeId, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, result);
            if (code < 0) {
                throw new GATKException(String.format("getting doubles from data-set '%s' in file '%s' resulted in code: %d", fullPath, file, code));
            }
            return result;
        });
    }

    /**
     * Reads a double matrix from a particular position in the underlying HDF5 file.
     *
     * @param fullPath the path.
     * @return never {@code null}, non-existing data or the wrong type in the HDF5
     * file will result in a {@link GATKException} instead.
     * @throws IllegalArgumentException if {@code fullPath} is {@code null}.
     * @throws GATKException if {@code fullPath} does not exist, contains the wrong data type (non-double) or
     *    its dimension is not 2.
     */
    public double[][] readDoubleMatrix(final String fullPath) {
        return readDataset(fullPath, (dataSetId, typeId, dimensions) -> {
            if (dimensions.length != 2) {
                throw new GATKException(
                        String.format("expected 2D double matrix for data-set '%s' in '%s' but it is %d-D", fullPath, file, dimensions.length));
            }
            final int rows = (int) dimensions[0];
            final int columns = (int) dimensions[1];
            final int size = rows * columns;
            final double[] values = new double[size];
            final int code = H5.H5Dread_double(dataSetId, typeId, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, values);
            if (code < 0) {
                throw new GATKException(String.format("getting double matrix from data-set '%s' in file '%s' resulted in code: %d", fullPath, file, code));
            }
            final double[][] result = new double[rows][columns];
            for (int i = 0; i < result.length; i++) {
                System.arraycopy(values, i * columns, result[i], 0, columns);
            }
            return result;
        });
    }

    /**
     * Functional interface for lambda used to read a HDF5 data-set using {@link #readDataset}.
     *
     * @param <T>
     */
    @FunctionalInterface
    private interface DatasetReader<T> {
        /**
         * Reads and returns the data in the underlying HDF5 file given its data-set id, type-id and dimensions.
         * @param dataSetId the target data-set.
         * @param typeId the target data-set's type.
         * @param dimensions the dimensions of the data-set.
         * @return might be {@code null}, is up to the implementation.
         * @throws HDF5LibraryException forwarding errors originated in the HD5F library.
         * @throws GATKException for any exceptional circumstance that make the returned value invalid (e.g. unexpected
         *  data-type or dimensions).
         */
        T apply(int dataSetId, int typeId, long[] dimensions) throws HDF5LibraryException;
    }

    /**
     * Dataset read template.
     *
     * <p>
     *     Reads the dataset located at the path provided, using the input datasetReader lambda.
     * </p>
     *
     * <p>
     *     This method takes care of allocating HDF5 data structures (Dataset, DataType and DataSpace)
     *     and freeing them after the data has been read by the lambda provided.
     * </p>
     *
     * <p>
     *     This method will attempt the closure of resources in the event that any exception occurred within
     *     the provided lambda.
     * </p>
     *
     * <p>
     *     The data-reading lambda may throw exceptions such as {@link GATKException} or {@link HDF5LibraryException}
     *     to indicate issues with
     *     the data (wrong type, wrong dimensions, etc.).
     * </p>
     *
     * @param fullPath the dataset full path.
     * @param datasetReader the actual reading operation implementation.
     * @param <T> the return data-type.
     * @return whatever {@code datasetReader} provided may return.
     * @throws IllegalArgumentException if either {@code fullPath} or {@code datasetReader} is {@code null}.
     * @throws GATKException if there is any error when reading the data, for example exceptions thrown by {@code datasetReader} or
     *    the HDF5 library.
     */
    private <T> T readDataset(final String fullPath, final DatasetReader<T> datasetReader) {
        if (fullPath == null) {
            throw new IllegalArgumentException("the path cannot be null");
        }
        checkIsOpen();
        int typeId = -1;
        int dataSetId = -1;
        int dataSpaceId = -1;
        try {
            dataSetId = openDataset(fullPath);
            typeId = openType(fullPath, dataSetId);
            dataSpaceId = openDataSpace(fullPath, dataSetId);
            final int dimNum = H5.H5Sget_simple_extent_ndims(dataSpaceId);
            final long[] dimensions = new long[dimNum];
            H5.H5Sget_simple_extent_dims(dataSpaceId, dimensions, null);
            return datasetReader.apply(dataSetId,typeId,dimensions);
        } catch (final HDF5LibraryException ex) {
            throw new GATKException(String.format("exception when reading from data-set '%s' in file '%s': %s", fullPath, file, ex.getMessage()), ex);
        } finally {
            closeResources(fullPath, typeId, dataSetId, dataSpaceId);
        }
    }

    /**
     * Closes HDF5 resources open during reading.
     *
     * @param fullPath the data-set full-path used in error messages.
     * @param typeId the data-type id.
     * @param dataSetId the data-set id.
     * @param dataSpaceId the data-space id.
     *
     * @throws GATKException if some exception occurred when closing the resources. If multiple errors took place,
     *    only one is reported as the cause of this exception.
     * @throws Error immediately if any was thrown during closure of any resource.
     */
    private void closeResources(final String fullPath, final int typeId, final int dataSetId, final int dataSpaceId) {

        final Optional<Exception> closureError = Stream.of(
                closeResource(() -> H5.H5Tclose(typeId)),
                closeResource(() -> H5.H5Sclose(dataSpaceId)),
                closeResource(() -> H5.H5Dclose(dataSetId)))
                .filter(Objects::nonNull).findFirst();

        if (closureError.isPresent()) {
            throw new GATKException(
                    String.format("exception when closing string retrieval on data-set '%s' in file '%s'", fullPath, file), closureError.get());
        }
    }

    /**
     * Common functional interface for resource closure actions.
     */
    @FunctionalInterface
    private interface ClosureAction {

         void run() throws Exception;
    }

    /**
     * Closes a HDF5 resource and captures exception thrown by the underlying HDF5 library.
     *
     * @param closure action require to close the resource.
     * @return reference to the exception if any was thrown, otherwise {@code null}.
     * @throws Error if any occurred.
     */
    private Exception closeResource(final ClosureAction closure) {
        try {
            closure.run();
        } catch (final Exception ex) {
            return ex;
        }
        return null;

    }

    /**
     * Opens a HDF5 DataSpace.
     *
     * @param fullPath dataset full path.
     * @param dataSetId dataset id.
     * @return the dataSpace id.
     * @throws HDF5LibraryException if there was some error thrown by the HDF5 library.
     * @throws GATKException if the HDF5 library returned a invalid dataSpace id indicating some kind of issue.
     */
    private int openDataSpace(final String fullPath, final int dataSetId) throws HDF5LibraryException {
        final int dataSpaceId = H5.H5Dget_space(dataSetId);
        if (dataSpaceId <= 0) {
            throw new GATKException(
                    String.format("getting the data-space of data-set '%s' in file '%s' resulted in code: %d", fullPath, file, dataSpaceId));
        }
        return dataSpaceId;
    }

    /**
     * Opens a HDF5 dataType.
     *
     * @param fullPath dataset full path.
     * @param dataSetId dataset id.
     * @return the dataType id.
     * @throws HDF5LibraryException if there was some error thrown by the HDF5 library.
     * @throws GATKException if the HDF5 library returned a invalid dataType id indicating some kind of issue.
     */
    private int openType(final String fullPath, final int dataSetId) throws HDF5LibraryException {
        final int typeId = H5.H5Dget_type(dataSetId);
        if (typeId <= 0) {
            throw new GATKException(
                    String.format("getting the type of data-set '%s' in file '%s' resulted in code: %d", fullPath, file, typeId));
        }
        return typeId;
    }

    /**
     * Opens a HDF5 dataSet.
     *
     * @param fullPath dataset full path.
     * @return the dataSet id.
     * @throws HDF5LibraryException if there was some error thrown by the HDF5 library.
     * @throws GATKException if the HDF5 library returned a invalid dataSet id indicating some kind of issue.
     */
    private int openDataset(final String fullPath) throws HDF5LibraryException {
        final int dataSetId = H5.H5Dopen(fileId, fullPath);
        if (dataSetId <= 0) {
            throw new GATKException(
                    String.format("opening string data-set '%s' in file '%s' failed with code: %d", fullPath, file, dataSetId));
        }
        return dataSetId;
    }

    /**
     * Checks that the reader is still open.
     *
     * <p>
     *     This check should be performed at the beginning of any reading operation.
     * </p>
     *
     * @throws IllegalStateException if the reader is closed.
     */
    private void checkIsOpen() {
        if (isClosed()) {
            throw new IllegalStateException("the reader is already closed");
        }
    }
}
