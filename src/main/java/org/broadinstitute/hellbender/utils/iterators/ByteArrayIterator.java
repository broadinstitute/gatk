package org.broadinstitute.hellbender.utils.iterators;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Trivial adapter class allowing a primitive byte[] array to be accessed using the java.util.Iterator interface
 */
public final class ByteArrayIterator implements Iterator<Byte> {
    private final byte[] byteArray;
    private int currentPosition;
    private int stopIndex;

    public ByteArrayIterator( final byte[] byteArray ) {
        this(byteArray, 0, byteArray.length);
    }

    /**
     * Iterates through this array, starting at index 'firstIndex' and going until
     * (but not including) index 'stopIndex'.
     */
    public ByteArrayIterator( final byte[] byteArray, int firstIndex, int stopIndex ) {
        if (!(byteArray.length==0 && firstIndex==0)) {
            // special case: empty buffer, we can start at 0 even though it's technically outside.
            Utils.validIndex(firstIndex, byteArray.length);
        }
        Utils.validateArg(stopIndex <= byteArray.length, () -> "stopIndex is "+stopIndex+" yet we only have "+byteArray.length+" bytes.");
        Utils.validateArg(stopIndex >= firstIndex, () -> "stopIndex<firstIndex ("+(stopIndex)+"<"+firstIndex+")");
        this.byteArray = byteArray;
        currentPosition = firstIndex;
        this.stopIndex = stopIndex;
    }

    @Override
    public boolean hasNext() {
        return currentPosition < stopIndex;
    }

    @Override
    public Byte next() {
        if (currentPosition >= stopIndex) {
            throw new NoSuchElementException("No more elements in byte array");
        }

        return byteArray[currentPosition++];
    }
}
