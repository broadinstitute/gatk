package org.broadinstitute.hellbender.utils.iterators;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Trivial adapter class allowing a primitive byte[] array to be accessed using the java.util.Iterator interface
 */
public final class ByteArrayIterator implements Iterator<Byte> {
    private final byte[] byteArray;
    private int currentPosition;

    public ByteArrayIterator( final byte[] byteArray ) {
        this.byteArray = byteArray;
        currentPosition = 0;
    }

    public boolean hasNext() {
        return currentPosition < byteArray.length;
    }

    public Byte next() {
        if ( currentPosition >= byteArray.length ) {
            throw new NoSuchElementException("No more elements in byte array");
        }

        return byteArray[currentPosition++];
    }
}
