/**
 * $Id: NybbleMatrix.java 72512 2008-09-12 17:16:43Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils.align;

/**
 * A two-dimensional matrix of nybbles.  Boundaries are checked only by assertions,
 * so don't count on run-time exceptions to save you if you call get or set with
 * bogus indices.
 *
 * @author tsharpe
 * @version $Revision$
 */
public final class NybbleMatrix
{
    public NybbleMatrix( int nCols, int nRows )
    {
        mNCols = nCols;
        mNRows = nRows;
        mVals = new byte[(nCols*nRows+1)/2];
    }

    /**
     * Get the value of the nybble.
     * @param col 0 <= col < nCols
     * @param row 0 <= row < nRows
     * @return The nybble.
     */
    public int get( int col, int row )
    {
        assert( col >= 0 && col < mNCols && row >= 0 && row < mNRows );

        int idx = row * mNCols + col;
        int val = mVals[idx/2];
        if ( (idx & 1) == 1 )
        {
            val >>= 4;
        }
        return val & 0x0F;
    }

    /**
     * Set the value of the nybble.
     * @param col 0 <= col < nCols
     * @param row 0 <= row < nRows
     * @param val 0 <= val < 16
     */
    public void set( int col, int row, int val )
    {
        assert( col >= 0 && col < mNCols && row >= 0 && row < mNRows );
        assert( val >= 0 && val < 16 );

        int idx = row * mNCols + col;
        int offset = idx / 2;
        int oldVal = mVals[offset];
        if ( (idx & 1) == 1 )
        {
            mVals[offset] = (byte)((oldVal & 0x0F) | (val << 4));
        }
        else
        {
            mVals[offset] = (byte)((oldVal & 0xF0) | (val & 0x0F));
        }
    }

    private int mNCols;
    private int mNRows;
    byte[] mVals;

    public final class SerialFiller
    {
        public void set( int val )
        {
            if ( mHigh )
            {
                mHigh = false;
                mVals[mIdx++] = (byte)((val << 4) | (mVal & 0x0F));
            }
            else
            {
                mHigh = true;
                mVal = val;
            }
        }

        public void done()
        {
            if ( mHigh )
            {
                mVals[mIdx] = (byte)(mVal & 0x0F);
            }
        }

        private int mIdx;
        private int mVal;
        private boolean mHigh;
    }
}
