/**
 * $Id: Score.java 57447 2008-01-28 17:41:24Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils.align;

/**
 * The score of an optimal alignment.
 * Also tracks the ending position of the best score w.r.t. the two sequences aligned.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class Score
{
    Score()
    {
        mScore = Float.NEGATIVE_INFINITY;
        mXBestIdx = -1;
        mYBestIdx = -1;
    }

    Score( Score score )
    {
        mScore = score.mScore;
        mXBestIdx = score.mXBestIdx;
        mYBestIdx = score.mYBestIdx;
    }

    /**
     * The best (largest) score.
     */
    public float getScore()
    {
        return mScore;
    }

    /**
     * The X position of the best (largest) score.
     */
    public int getXBestIndex()
    {
        return mXBestIdx;
    }

    /**
     * The Y position of the best (largest) score.
     */
    public int getYBestIndex()
    {
        return mYBestIdx;
    }

    protected void checkRowEnd( float score, int xIdx, int yIdx )
    {
        checkScore(score,xIdx,yIdx);
    }

    protected void checkLastRow( float[] scores, int yIdx )
    {
        for ( int xIdx = 0; xIdx < scores.length; ++xIdx )
        {
            checkScore(scores[xIdx],xIdx,yIdx);
        }
    }

    protected final void checkScore( float score, int xIdx, int yIdx )
    {
        if ( score > mScore )
        {
            mScore = score;
            mXBestIdx = xIdx;
            mYBestIdx = yIdx;
        }
    }

    private float mScore;
    private int mXBestIdx;
    private int mYBestIdx;
}
