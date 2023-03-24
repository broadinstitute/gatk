/*
 * $Id: SingletonIterator.java 49074 2007-10-02 20:06:53Z tsharpe $
 * WHITEHEAD INSTITUTE
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2002 by the
 * Whitehead Institute for Biomedical Research.  All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support
 * whatsoever.  The Whitehead Institute can not be responsible for its
 * use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * An iterator over a single object.
 *
 * @author tsharpe
 * @version $Revision: 28397 $
 */
public class SingletonIterator<T>
    implements Iterator<T>
{
    /**
     * Make one.  If obj is null, then hasNext will return false immediately.
     * @param obj  The object you'd like the iterator to present.
     */
    SingletonIterator( T obj )
    {
        mObj = obj;
    }

    /**
     * Not supported.
     */
    public void remove()
    {
        throw new UnsupportedOperationException("Read-only iterator.");
    }

    public boolean hasNext()
    {
        return mObj != null;
    }

    public T next()
    {
        T result = mObj;
        if ( result == null )
        {
            throw new NoSuchElementException("Iteration complete.");
        }
        mObj = null;
        return result;
    }

    private T mObj;
}
