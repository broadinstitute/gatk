package org.broadinstitute.hellbender.utils.collections;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class NestedIntegerArray<T> {

    private static Logger logger = LogManager.getLogger(NestedIntegerArray.class);

    protected final Object[] data;

    protected final int numDimensions;
    protected final int[] dimensions;

    // Preallocate the first two dimensions to limit contention during tree traversals in put()
    private static final int NUM_DIMENSIONS_TO_PREALLOCATE = 2;

    public NestedIntegerArray(final int... dimensions) {
        numDimensions = dimensions.length;
        if ( numDimensions == 0 )
            throw new IllegalArgumentException("There must be at least one dimension to an NestedIntegerArray");
        this.dimensions = Arrays.copyOf(dimensions, dimensions.length);

        int dimensionsToPreallocate = Math.min(dimensions.length, NUM_DIMENSIONS_TO_PREALLOCATE);

        if ( logger.isDebugEnabled() ) logger.debug(String.format("Creating NestedIntegerArray with dimensions %s", Arrays.toString(dimensions)));
        if ( logger.isDebugEnabled() ) logger.debug(String.format("Pre-allocating first %d dimensions", dimensionsToPreallocate));

        data = new Object[dimensions[0]];
        preallocateArray(data, 0, dimensionsToPreallocate);

        if ( logger.isDebugEnabled() ) logger.debug(String.format("Done pre-allocating first %d dimensions", dimensionsToPreallocate));
    }

    /**
     * @return the dimensions of this nested integer array.  DO NOT MODIFY
     */
    public int[] getDimensions() {
        return dimensions;
    }

    /**
     * Recursively allocate the first dimensionsToPreallocate dimensions of the tree
     *
     * Pre-allocating the first few dimensions helps limit contention during tree traversals in put()
     *
     * @param subarray current node in the tree
     * @param dimension current level in the tree
     * @param dimensionsToPreallocate preallocate only this many dimensions (starting from the first)
     */
    private void preallocateArray( Object[] subarray, int dimension, int dimensionsToPreallocate ) {
        if ( dimension >= dimensionsToPreallocate - 1 ) {
            return;
        }

        for ( int i = 0; i < subarray.length; i++ ) {
            subarray[i] = new Object[dimensions[dimension + 1]];
            preallocateArray((Object[])subarray[i], dimension + 1, dimensionsToPreallocate);
        }
    }

    @SuppressWarnings("unchecked")
    public T get(final int... keys) {
        final int numNestedDimensions = numDimensions - 1;
        Object[] myData = data;

        for( int i = 0; i < numNestedDimensions; i++ ) {
            if ( keys[i] >= dimensions[i] )
                return null;

            myData = (Object[])myData[keys[i]];
            if ( myData == null )
                return null;
        }

        return (T)myData[keys[numNestedDimensions]];
    }

    /**
     * Insert a value at the position specified by the given keys.
     *
     * @param value value to insert
     * @param keys keys specifying the location of the value in the tree
     */
    public void put(final T value, final int... keys) { // WARNING! value comes before the keys!
        if ( keys.length != numDimensions )
            throw new IllegalArgumentException("Exactly " + numDimensions + " keys should be passed to this NestedIntegerArray but " + keys.length + " were provided");

        final int numNestedDimensions = numDimensions - 1;
        Object[] myData = data;
        for ( int i = 0; i < numNestedDimensions; i++ ) {
            if ( keys[i] >= dimensions[i] )
                throw new IllegalArgumentException("Key " + keys[i] + " is too large for dimension " + i + " (max is " + (dimensions[i]-1) + ")");

            // If we're at or beyond the last dimension that was pre-allocated, we need to
            // check to see if the next branch exists, and if it doesn't, create it
            if ( i >= NUM_DIMENSIONS_TO_PREALLOCATE - 1 ) {
                if ( myData[keys[i]] == null ) {
                    myData[keys[i]] = new Object[dimensions[i + 1]];
                }
            }

            myData = (Object[])myData[keys[i]];
        }

        myData[keys[numNestedDimensions]] = value;
    }

    public List<T> getAllValues() {
        final List<T> result = new ArrayList<>();
        fillAllValues(data, result);
        return result;
    }

    @SuppressWarnings("unchecked")
    private void fillAllValues(final Object[] array, final List<T> result) {
        for ( Object value : array ) {
            if ( value == null )
                continue;
            if ( value instanceof Object[] )
                fillAllValues((Object[])value, result);
            else
                result.add((T)value);
        }
    }

    public static class Leaf<T> {
        public final int[] keys;
        public final T value;

        public Leaf(final int[] keys, final T value) {
            this.keys = keys;
            this.value = value;
        }
    }

    public List<Leaf<T>> getAllLeaves() {
        final List<Leaf<T>> result = new ArrayList<>();
        fillAllLeaves(data, new int[0], result);
        return result;
    }

    @SuppressWarnings("unchecked")
    private void fillAllLeaves(final Object[] array, final int[] path, final List<Leaf<T>> result) {
        for ( int key = 0; key < array.length; key++ ) {
            final Object value = array[key];
            if ( value == null )
                continue;
            final int[] newPath = appendToPath(path, key);
            if ( value instanceof Object[] ) {
                fillAllLeaves((Object[]) value, newPath, result);
            } else {
                result.add(new Leaf<T>(newPath, (T)value));
            }
        }
    }

    private int[] appendToPath(final int[] path, final int newKey) {
        final int[] newPath = new int[path.length + 1];
        for ( int i = 0; i < path.length; i++ )
            newPath[i] = path[i];
        newPath[path.length] = newKey;
        return newPath;
    }
}
