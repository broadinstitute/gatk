/*
 * $Id: DoubleArrayList.java 70611 2008-07-29 18:06:55Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;


/**
 * A series of doubles supported by an array.
 *
 * @author tsharpe
 * @version $Revision: 70611 $
 */
public class DoubleArrayList
    implements NumberList
{
    public DoubleArrayList()
    {
        this(DEFAULT_INITIAL_CAPACITY);
    }

    public DoubleArrayList( int capacity )
    {
        mVals = new double[capacity];
    }

    public DoubleArrayList( double[] vals )
    {
        mVals = vals == null ? NULL_LIST : vals;
        mSize = mVals.length;
    }

    public int size()
    {
        return mSize;
    }

    public int intVal( int idx )
    {
        return Math.round((float)dblVal(idx));
    }

    public double dblVal( int idx )
    {
        if ( idx < 0 )
        {
            throw new IndexOutOfBoundsException(idx + "<0");
        }
        if ( idx >= mSize )
        {
            throw new IndexOutOfBoundsException(idx + ">=" + mSize);
        }
        return mVals[idx];
    }

    public synchronized void add( double val )
    {
        if ( mSize >= mVals.length )
        {
            double[] newVals = new double[mVals.length*2 + 5];
            System.arraycopy(mVals,0,newVals,0,mVals.length);
            mVals = newVals;
        }
        mVals[mSize++] = val;
    }

    public Iterator iterator()
    {
        return new NumberListUtil.Itr(this);
    }

    private double[] mVals;
    private int mSize;
    private static final double[] NULL_LIST = new double[0];
    private static final int DEFAULT_INITIAL_CAPACITY = 15;
}
