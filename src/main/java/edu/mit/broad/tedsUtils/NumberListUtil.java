/**
 * $Id: NumberListUtil.java 72259 2008-09-03 20:39:10Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils;

import java.util.NoSuchElementException;

import edu.mit.broad.tedsUtils.NumberList.Iterator;


/**
 * Utilities for NumberLists.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class NumberListUtil
{
    public static class SubList
        implements NumberList
    {
        public SubList( NumberList list, int start, int end )
        {
            if ( list == null )
            {
                throw new IllegalArgumentException("NumberList cannot be null.");
            }
            if ( start < 0 || start > list.size() )
            {
                throw new IllegalArgumentException("Start out of bounds.");
            }
            if ( end < start || end > list.size() )
            {
                throw new IllegalArgumentException("End out of bounds.");
            }
            mList = list;
            mStart = start;
            mEnd = end;
        }

        public double dblVal( int idx )
        {
            if ( idx < 0 || idx > size() )
            {
                throw new IllegalArgumentException("Index is out of bounds.");
            }
            return mList.dblVal(mStart+idx);
        }

        public int intVal( int idx )
        {
            if ( idx < 0 || idx > size() )
            {
                throw new IllegalArgumentException("Index is out of bounds.");
            }
            return mList.intVal(idx);
        }

        public Iterator iterator()
        {
            return new Itr(this);
        }

        public int size()
        {
            return mEnd - mStart;
        }

        private NumberList mList;
        private int mStart;
        private int mEnd;
    }

    public static class Itr
        implements Iterator
    {
        public Itr( NumberList list )
        {
            mList = list;
        }

        public boolean hasNext()
        {
            return mIdx < mList.size();
        }

        public int nextInt()
        {
            if ( mIdx > mList.size() )
            {
                throw new NoSuchElementException("We're off the end of the iterator.");
            }
            return mList.intVal(mIdx++);
        }

        public double nextDbl()
        {
            if ( mIdx > mList.size() )
            {
                throw new NoSuchElementException("We're off the end of the iterator.");
            }
            return mList.dblVal(mIdx++);
        }

        private NumberList mList;
        private int mIdx;
    }

    public static class MovingAverage
    {
        public MovingAverage( NumberList list, int start, int length )
        {
            if ( list == null )
            {
                throw new IllegalArgumentException("number list cannot be null");
            }
            if ( start < 0 || length < 1 || start + length > list.size() )
            {
                throw new IllegalArgumentException("bad indexes");
            }

            mList = list;
            mStart = start;
            mEnd = mStart + length;

            mTotal = 0.;
            for ( int idx = mStart; idx < mEnd; ++idx )
            {
                mTotal += mList.dblVal(idx);
            }
        }

        public double total()
        {
            return mTotal;
        }

        public double average()
        {
            return mTotal / (mEnd - mStart);
        }

        public int start()
        {
            return mStart;
        }

        public int end()
        {
            return mEnd;
        }

        public int center()
        {
            return (mStart + mEnd)/2;
        }

        public int length()
        {
            return mEnd - mStart;
        }

        public boolean next()
        {
            boolean result = false;
            if ( mEnd < mList.size() )
            {
                mTotal -= mList.dblVal(mStart++);
                mTotal += mList.dblVal(mEnd++);
                result = true;
            }
            return result;
        }

        public boolean prev()
        {
            boolean result = false;
            if ( mStart > 0 )
            {
                mTotal -= mList.dblVal(--mEnd);
                mTotal += mList.dblVal(--mStart);
                result = true;
            }
            return result;
        }

        public boolean incrementStart()
        {
            boolean result = false;
            if ( mStart + 1 < mEnd )
            {
                mTotal -= mList.dblVal(mStart++);
                result = true;
            }
            return result;
        }

        public boolean decrementStart()
        {
            boolean result = false;
            if ( mStart > 0 )
            {
                mTotal += mList.dblVal(--mStart);
                result = true;
            }
            return result;
        }

        public boolean incrementEnd()
        {
            boolean result = false;
            if ( mEnd < mList.size() )
            {
                mTotal += mList.dblVal(mEnd++);
                result = true;
            }
            return result;
        }

        public boolean decrementEnd()
        {
            boolean result = false;
            if ( mEnd - 1 > mStart )
            {
                mTotal -= mList.dblVal(--mEnd);
                result = true;
            }
            return result;
        }

        private NumberList mList;
        private int mStart;
        private int mEnd;
        private double mTotal;
    }

    public static double average( NumberList list )
    {
        int nnn = list.size();
        double tot = 0.;
        for ( int idx = 0; idx < nnn; ++idx )
        {
            tot += list.dblVal(idx);
        }
        return tot/nnn;
    }
    public static int maxArg( NumberList list )
    {
        double maxVal = -Double.MAX_VALUE;
        int result = 0;
        int nnn = list.size();
        for ( int idx = 0; idx < nnn; ++idx )
        {
            double val = list.dblVal(idx);
            if ( val > maxVal )
            {
                maxVal = val;
                result = idx;
            }
        }
        return result;
    }
}
