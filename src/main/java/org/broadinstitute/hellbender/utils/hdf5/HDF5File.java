package org.broadinstitute.hellbender.utils.hdf5;

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.HDF5Constants;
import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;
import ncsa.hdf.hdf5lib.exceptions.HDF5LibraryException;
import ncsa.hdf.hdf5lib.structs.H5G_info_t;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.IntSupplier;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * HDF5 File reader.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class HDF5File implements AutoCloseable {

    /**
     * Dimensions HDF5 object used for single scalar value data-sets.
     */
    private static final long[] SCALAR_VALUE_DIMENSIONS = new long[] { 1 };

    /**
     * String used to separated elements in HDF5 path.
     */
    private static final String PATH_ELEMENT_SEPARATOR = "/";

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
     * Indicates whether this file can be written.
     */
    private final boolean canWrite;

    /**
     * Special {@link #fileId} value to indicate that the reader is closed.
     */
    protected final int FILE_ID_WHEN_CLOSED = -1;

    /**
     * Creates a new HDF5-file accessor object given the file name for reading purposes.
     * @param file the underlying file name.
     */
    public HDF5File(final File file) {
        this(file, OpenMode.READ_ONLY);
    }

    /**
     * Creates a new HDF5 reader on a existing file.
     *
     * @param file the target file.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws GATKException if the HDF5 library is not supported or could not be initialized.
     */
    public HDF5File(final File file, final OpenMode mode) {
        HDF5Library.getLibrary();
        fileId = open(this.file = Utils.nonNull(file, "the input file cannot be null"), Utils.nonNull(mode, "the mode cannot be null"));
        if (fileId < 0) {
            throw new GATKException(
                    String.format("failure when opening '%s' for read-only access; negative fileId: %d",file.getAbsolutePath(),fileId)
            );
        }
        canWrite = mode.canWrite();
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
            final boolean isVariableString = H5.H5Tis_variable_str(typeId);

            final String[] result = new String[(int) dimensions[0]];

            final int code;
            if (isVariableString) {
                code = H5.H5DreadVL(dataSetId, typeId, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, result);
            }
            else {
                code = H5.H5Dread_string(dataSetId, typeId, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, result);
            }
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
        final int dataSetId = H5.H5Dopen(fileId, fullPath, HDF5Constants.H5P_DEFAULT);
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

    /**
     * Opens a HDF5 file given its access {@link OpenMode mode}.
     */
    private static int open(final File file, final OpenMode mode) {
        final int fileId;
        try {
            if (mode == OpenMode.CREATE) {
                file.delete();
                file.createNewFile();
            }
            fileId = H5.H5Fopen(file.getAbsolutePath(), mode.flags, HDF5Constants.H5P_DEFAULT);
        } catch (final HDF5LibraryException | IOException e) {
            throw new GATKException(
                    String.format("exception when opening '%s' with %s mode: %s",file.getAbsolutePath(), mode, e.getMessage()), e);
        }
        if (fileId < 0) {
            throw new GATKException(
                    String.format("failure when opening '%s' for read-only access; negative fileId: %d",file.getAbsolutePath(),fileId)
            );
        }
        return fileId;
    }

    /**
     * Different ways a HDF5 file can be opened.
     */
    public enum OpenMode {
        /**
         * Access the file for read-only, it won't try to create an empty file if none exists.
         */
        READ_ONLY(HDF5Constants.H5F_ACC_RDONLY),
        /**
         * Access the file for read or write, it won't try to create an empty file if non exists.
         */
        READ_WRITE(HDF5Constants.H5F_ACC_RDWR),

        /**
         * Creates a new file with empty contents and give read-write access to it.
         * <p>
         * It will overwrite its contents.
         * </p>
         */
        CREATE(HDF5Constants.H5F_ACC_RDWR);

        private final int flags;

        OpenMode(final int flags) {
            this.flags = flags;
        }

        /**
         * Checks whether this mode allows for writing operations.
         * @return {@code true} iff the mode allows for writing operations.
         */
        public boolean canWrite() {
            return READ_ONLY != this;
        }
    }

    //////////////////////
    // Write operations.
    //////////////////////

    /**
     * Create dataset group (directory) in the HD5 file.
     * <p>
     * This method ensures that the requested group exist after this invocation creating
     * any required parent group recursively.
     * </p>
     * <p>
     * Existing groups won't be modified.
     * </p>
     * <p>
     * The path elements are separated by the slash character.
     * </p>
     * @param path absolute path for the group to create.
     * @throws IllegalArgumentException if {@code path} is {@code null}.
     * @throws GATKException if there was some issue trying to create the group. This includes but is not limited to:
     *     some low level exception if the the library, an existent non-group object standing in
     * @return {@code true} iff this operation actually modified the file, i.e. the group or any of its parent didn't exists.
     */
    public boolean makeGroup(final String path) throws GATKException {
        Utils.nonNull(path, "the path cannot be null");
        checkCanWrite();
        boolean modified = false;
        final Queue<String> pathElements = Stream.of(path.split(PATH_ELEMENT_SEPARATOR)).filter(s -> !s.isEmpty()).collect(Collectors.toCollection(ArrayDeque::new));
        final ArrayList<Integer> groupIds = new ArrayList<>(pathElements.size() + 1);
        final List<String> pathSoFar = new ArrayList<>(pathElements.size());
        try {
            final int rootId = H5.H5Gopen(fileId, PATH_ELEMENT_SEPARATOR, HDF5Constants.H5P_DEFAULT);
            if (rootId < 0) {
                throw new GATKException(String.format("there was a problem to find a group (%s) in file %s", path, file));
            }

            groupIds.add(rootId);
            while (!pathElements.isEmpty()) {
                final int parentId = groupIds.get(groupIds.size() - 1);
                final String childName = pathElements.remove();
                final int childType = findOutGroupChildType(parentId, childName, path);
                if (childType == HDF5Constants.H5G_UNKNOWN) {
                    modified = true;
                    int nextGroupId = H5.H5Gcreate(parentId, childName, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT);
                    if (nextGroupId < 0) {
                        throw new GATKException(String.format("there was a problem to make group (%s) in file %s", path, file));
                    }
                    groupIds.add(nextGroupId);
                } else if (childType != HDF5Constants.H5G_GROUP) {
                    throw new GATKException(String.format("there is a non-group object with type (%d) on the way to the requested group name (%s) in file %s", childType, path, file));
                } else {
                    final int nextGroupId = H5.H5Gopen(parentId, childName, HDF5Constants.H5P_DEFAULT);
                    if (nextGroupId < 0) {
                        throw new GATKException(String.format("problem trying to find a group (%s) in file %s", path, file));
                    }
                    groupIds.add(nextGroupId);
                }
                pathSoFar.add(childName);
            }

        } catch (final HDF5LibraryException ex) {
            throw new GATKException(String.format("Exception trying to find/make a group (%s) in file %s", path, file), ex);
        } finally {
            // we need to do our best to close all the groups we have opened.
            for (final int groupId : groupIds) {
                try { H5.H5Gclose(groupId); } catch (final HDF5LibraryException ex) { }
            }
        }
        return modified;
    }

    /**
     * Returns the type of a group child given.
     * <p>
     * Type constants are listed in {@link HDF5Constants}, eg.:
     * <ul>
     *     <li>{@link HDF5Constants#H5G_GROUP H5G_GROUP} indicate that the child is another group</li>
     *     <li>{@link HDF5Constants#H5G_DATASET H5G_DATASET} indicate that the child is a dataset...</li>
     * </ul>
     *
     * </p>
     * {@link HDF5Constants#H5G_UNKNOWN H5G_UNKNOWN} indicates that the child is not present.
     *
     * @param groupId the parent group id. It must be open.
     * @param name of the target child.
     * @param fullPath full path reported in exceptions when there is an issue.
     * @return {@link HDF5Constants#H5G_UNKNOWN H5G_UNKNOWN} if there is no such a child node, other wise any other valid type constant.
     * @throws HDF5LibraryException if any is thrown by the HDF5 library.
     */
    private int findOutGroupChildType(final int groupId, final String name, final String fullPath) throws HDF5LibraryException {

        // Use an single position array to return a value is kinda inefficient but that is the way it is:
        final long[] numObjsResult = new long[1];
        H5G_info_t result = H5.H5Gget_info(groupId);
        numObjsResult[0] = result.nlinks;

        final int childCount = (int) numObjsResult[0];
        if (childCount == 0) { // this is no premature optimization: get_obj_info_all really cannot handle length 0 arrays.
            return HDF5Constants.H5G_UNKNOWN;
        } else {
            final String[] childNames = new String[childCount];
            final int[] childTypes = new int[childCount];
            final int[] lTypes = new int[childCount];
            final long[] childRefs = new long[childCount];

            // Example call in HDF docs (https://www.hdfgroup.org/HDF5/examples/api18-java.html .... H5_Ex_G_Iterate.java Line 71):
            //  H5.H5Gget_obj_info_all(file_id, DATASETNAME, oname, otype, ltype, orefs, HDF5Constants.H5_INDEX_NAME);
            if (H5.H5Gget_obj_info_all(groupId, ".", childNames, childTypes, lTypes, childRefs, HDF5Constants.H5_INDEX_NAME) < 0) {
                throw new GATKException(String.format("problem trying to find a group (%s) in file %s", fullPath, file));
            }
            final int childIndex = ArrayUtils.indexOf(childNames, name);
            if (childIndex == -1) {
                return HDF5Constants.H5G_UNKNOWN;
            } else {
                return childTypes[childIndex];
            }
        }
    }

    /**
     * Creates or overwrites a double value at particular position in the underlying HDF5 file.
     *
     * @param fullPath the path where to place the double value in the HDF5 file.
     * @return true iff the new data-set had to be created (none existed for that path).
     * @throws IllegalArgumentException if {@code fullPath} is {@code null} or is not a valid data-type name.
     * @throws GATKException if {@code fullPath} does not exist, contains the wrong data type (non-double) or
     *    is an array with several values or multidimensional.
     */
    public boolean makeDouble(final String fullPath, final double value) {
        return makeDataset(fullPath, basicTypeCopyIdSupplier(HDF5Constants.H5T_INTEL_F64), SCALAR_VALUE_DIMENSIONS, new double[] {value});
    }

    /**
     * Creates or overwrites a double array at particular position in the underlying HDF5 file.
     *
     * @param fullPath the path where to place the double array in the HDF5 file.
     * @return true iff the new data-set had to be created (none existed for that path).
     * @throws IllegalArgumentException if {@code fullPath} is {@code null} or is not a valid data-type name.
     * @throws GATKException if {@code fullPath} does not exist, contains the wrong data type (non-double) or
     *    is not a 1D array or is too small to contain the new value.
     */
    public boolean makeDoubleArray(final String fullPath, final double[] value) {
        Utils.nonNull(value);
        final long[] dimensions = new long[] { value.length };
        return makeDataset(fullPath, basicTypeCopyIdSupplier(HDF5Constants.H5T_INTEL_F64), dimensions,  value);
    }

    /**
     * Creates or overwrites a double 2D matrix at particular position in the underlying HDF5 file.
     *
     * @param fullPath the path where to place the double matrix in the HDF5 file.
     * @return the stored value.
     * @throws IllegalArgumentException if {@code fullPath} is {@code null} or is not a valid data-type name.
     * @throws GATKException if {@code fullPath} does not exist, contains the wrong data type (non-double) or
     *    is not a 2D array or is too small to contain the new value.
     */
    public boolean makeDoubleMatrix(final String fullPath, final double[][] value) {
        Utils.nonNull(value, "the value provided cannot be null");
        if (value.length == 0) {
            throw new IllegalArgumentException("the value provided must have some elements");
        }
        final int columnCount = Utils.nonNull(value[0], "the input value array cannot contain nulls: 0").length;
        if (columnCount == 0) {
            throw new IllegalArgumentException("the value provided must have some elements");
        }
        for (int i = 1; i < value.length; i++) {
            if (Utils.nonNull(value[i], "some row data is null: " + i).length != columnCount) {
                throw new IllegalArgumentException("some rows in the input value matrix has different number of elements");
            }
        }
        final long[] dimensions = new long[] { value.length, columnCount };
        return makeDataset(fullPath, basicTypeCopyIdSupplier(HDF5Constants.H5T_INTEL_F64), dimensions, value);
    }

    /**
     * Creates or overwrites a string array dataset given its full path within the HDF5 file and values.
     * @param fullPath the full path of the dataset inside the HDF5 file.
     * @param values the string array to write in the file.
     * @return {@code true} if a new data-set need to be created.
     * @throws IllegalArgumentException if {@code values} is {@code null} or contains {@code null}.
     * @throws GATKException if there was any low-level issue accessing the HDF5 file.
     */
    public boolean makeStringArray(final String fullPath, final String ... values) {
        Utils.nonNull(values);
        final long[] dimensions = new long[] { values.length };
        for (final String value : values) {
            Utils.nonNull(value);
        }

        return makeDataset(fullPath, basicStringArrayTypeIdSupplier(), dimensions, values);
    }

    /**
     * Check that the result of a {@link H5} utility class call does not indicate an error.
     * <p>
     * If the result is less than 0, this indicates an error in the execution of the generating
     * call and a {@link GATKException} is thrown.
     * </p>
     * <p>
     * If the result is 0 or greater we return that value to the caller.
     * </p>
     * @param result the previous H5 call result code to check.
     * @param exceptionMessageSupplier supplies a message for the exception.
     * @return sames as {@code result}.
     */
    private int checkH5Result(final int result, final Supplier<String> exceptionMessageSupplier) {
        if (result < 0) {
            throw new GATKException(exceptionMessageSupplier != null ? exceptionMessageSupplier.get(): "");
        } else {
            return result;
        }
    }

    /**
     * Creates a HDF5 type based on a simple copy of a HDF5 type class.
     * @param classId the class id to copy.
     * @return never {@code null}.
     */
    private IntSupplier basicTypeCopyIdSupplier(final int classId) {
        return () -> {
            try {
                return checkH5Result(H5.H5Tcopy(classId), () -> "");
            } catch (final HDF5LibraryException ex) {
                throw new GATKException("", ex);
            }
        };
    }

    /**
     * Creates a HDF5 type suitable for string arrays.
     */
    private IntSupplier basicStringArrayTypeIdSupplier() {
        return () -> {
            int result = -1;
            try {
              try {
                  result = checkH5Result(H5.H5Tcopy(HDF5Constants.H5T_C_S1), () -> "problem copying string type id");
                  checkH5Result(H5.H5Tset_size(result, HDF5Constants.H5T_VARIABLE), () -> "problem setting maximum size of string type to ");
                  return result;
              } catch (final HDF5LibraryException ex) {
                  throw new GATKException("", ex);
              }
            } catch (final GATKException ex) {
                // if we fail somewhere after creating the type, we need to close it if we can:
                if (result != -1) { try { H5.H5Tclose(result); } catch (final HDF5Exception ex2) {} }
                throw ex;
            }
        };
    }

    /**
     * General dataset making recipe.
     * @param fullPath the dataset full path.
     * @param typeIdSupplier type id supplier lambda.
     * @param dimensions array with the dimensions of the data.
     * @param data the data. It must be an array of the appropriate type given the type that is
     *             going to be returned by the {@code typeIdSupplier}.
     * @return true iff the data-set needed to be created (it did not existed previously). It will
     * return false if the data-set existed even if it was modified in the process.
     */
    private boolean makeDataset(final String fullPath, final IntSupplier typeIdSupplier, final long[] dimensions, final Object data) {
        checkCanWrite();
        int typeCopyId = -1;
        try {
            typeCopyId = typeIdSupplier.getAsInt();
            final Pair<String, String> pathAndName = splitPathInParentAndName(fullPath);
            final String groupPath = pathAndName.getLeft();
            final String dataSetName = pathAndName.getRight();
            makeGroup(groupPath);
            final int childType = findOutGroupChildType(groupPath, dataSetName, fullPath);
            if (childType == HDF5Constants.H5G_UNKNOWN) {
                createDataset(fullPath, typeCopyId, dimensions);
                writeDataset(fullPath, typeCopyId, data);
                return true;
            } else if (childType == HDF5Constants.H5G_DATASET) {
                writeDataset(fullPath, typeCopyId, data);
                return false;
            } else {
                throw new GATKException(String.format("problem trying to write dataset %s in file %s: there is a collision with a non-dataset object", fullPath, file));
            }
        } finally {
            if (typeCopyId != -1) { try { H5.H5Tclose(typeCopyId); } catch (final HDF5Exception ex ){} }
        }
    }

    /**
     * Returns the a group's child type given the group path and the child name.
     * @param groupPath the group's full path
     * @param name the child name.
     * @param fullPath combination of groupPath and name, used for exception messages.
     * @return {@link HDF5Constants#H5G_UNKNOWN} if not such a child exist. {@link HDF5Constants#H5G_DATASET} for
     * dataset, {@link HDF5Constants#H5G_GROUP} for groups, etc...
     */
    private int findOutGroupChildType(final String groupPath, final String name, final String fullPath) {
        int groupId = -1;
        try {
            groupId = H5.H5Gopen(fileId, groupPath, HDF5Constants.H5P_DEFAULT);
            final int childType = findOutGroupChildType(groupId, name, fullPath);
            return childType;
        } catch (final HDF5Exception ex) {
            throw new GATKException(String.format("problem when trying to resolve element %s type in file %s", fullPath, file));
        } finally {
            if (groupId != -1) {
                try { H5.H5Gclose(groupId); } catch (final HDF5Exception ex) { } }
        }
    }

    /**
     * General dataset creating HDF5 library recipe.
     * @param fullPath the target full-path within the file.
     * @param typeId the data type id.
     * @param dimensions the data dimensions.
     */
    private void createDataset(final String fullPath, final int typeId, final long[] dimensions) {
        int dataSpaceId = -1;
        int dataSetId = -1;
        try {
            dataSpaceId = checkH5Result(H5.H5Screate_simple(dimensions.length, dimensions, null),
                    () -> String.format("problem trying to create dataset %s in file %s", fullPath, file));
            dataSetId = checkH5Result(H5.H5Dcreate(fileId, fullPath, typeId, dataSpaceId,
                            HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT),
                    () -> String.format("problem trying to create dataset %s in file %s", fullPath, file));
        } catch (final HDF5Exception ex) {
            throw new GATKException(String.format("problem trying to create dataset %s in file %s", fullPath, file), ex);
        } finally {
            if (dataSetId != -1) { try { H5.H5Dclose(dataSetId); } catch (final HDF5LibraryException ex) {} }
            if (dataSpaceId != -1) { try { H5.H5Sclose(dataSpaceId); } catch (final HDF5LibraryException ex) {} }
        }
    }

    /**
     * General dataset writing HDF5 library recipe.
     * @param fullPath the target full-path within the file.
     * @param typeId the dataset type id.
     * @param data the data to write.
     */
    private void writeDataset(final String fullPath, final int typeId, final Object data) {
        int dataSetId = -1;
        int dataTypeId = -1;
        try {
            dataSetId = checkH5Result(H5.H5Dopen(fileId, fullPath, HDF5Constants.H5P_DEFAULT), () -> String.format("problem opening dataset %s in file %s: ", fullPath, file));
            dataTypeId = H5.H5Dget_type(dataSetId);
            if (!H5.H5Tequal(dataTypeId, typeId)) {
                throw new GATKException(String.format("problem writing new data to existing dataset %s: type is incompatible", fullPath));
            }
            checkH5Result(H5.H5Dwrite(dataSetId, typeId, HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT, data, false),
                    () -> String.format("error trying to write data-set %s in file %s", fullPath, file));
        } catch (final HDF5Exception ex) {
            throw new GATKException(String.format("problem writing dataset %s in file %s", fullPath, file), ex);
        } finally {
            if (dataTypeId != -1) try { H5.H5Tclose(dataTypeId); } catch (final HDF5LibraryException ex) {}
            if (dataSetId != -1) try { H5.H5Dclose(dataSetId);} catch (final HDF5LibraryException ex) {}
        }
    }

    /**
     * Decomposes a HDF5 path into the including group are data-set name.
     * @param fullPath the path to decompose.
     * @return never {@code null}, an immutable pair.
     * @throws IllegalArgumentException if {@code fullPath} is {@code null} or it finishes with a slash, indicating that the referenced object is a group.
     */
    private Pair<String,String> splitPathInParentAndName(final String fullPath) {
        final int lastSlashIndex = fullPath.lastIndexOf(PATH_ELEMENT_SEPARATOR);
        if (lastSlashIndex == -1) {
            return new ImmutablePair<>(PATH_ELEMENT_SEPARATOR,fullPath);
        } else if (lastSlashIndex == fullPath.length() - 1) {
            throw new IllegalArgumentException(String.format("the path provided make reference to a group name (directory) as finished with an slash: '%s'", fullPath));
        } else {
            return new ImmutablePair<>(fullPath.substring(0,lastSlashIndex),fullPath.substring(lastSlashIndex + 1));
        }
    }

    /**
     * Checks whether this HDF5 file handle supports writting operations.
     */
    private void checkCanWrite() {
        if (!canWrite) {
            throw new UnsupportedOperationException("this HDF5 file handle is not able to write");
        }
    }

}
