package org.broadinstitute.hellbender.utils.iterators;

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
        this.byteArray = byteArray;
        currentPosition = 0;
        stopIndex = byteArray.length;
    }

    /**
     * Iterates through this array, starting at index 'firstIndex' and going until
     * (but not including) index 'stopIndex'.
     */
    public ByteArrayIterator( final byte[] byteArray, int firstIndex, int stopIndex ) {
        if (stopIndex>byteArray.length) {
            throw new IllegalArgumentException("stopIndex should be at most byteArray.length ("+byteArray.length+"), but it is "+stopIndex);
        }
        this.byteArray = byteArray;
        currentPosition = firstIndex;
        this.stopIndex = stopIndex;
    }

    public boolean hasNext() {
        return currentPosition < stopIndex;
    }

    public Byte next() {
        if (currentPosition >= stopIndex) {
            throw new NoSuchElementException("No more elements in byte array");
        }

        return byteArray[currentPosition++];
    }
}
