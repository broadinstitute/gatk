/**
 * $Id: Aligner.java 73459 2008-10-10 16:10:03Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils.align;

import java.util.LinkedList;

import edu.mit.broad.tedsUtils.CharSubSequence;
import edu.mit.broad.tedsUtils.align.Alignment.Block;
import edu.mit.broad.tedsUtils.align.NybbleMatrix.SerialFiller;

/**
 * Base class for Alignment strategies.
 *
 * @author tsharpe
 * @version $Revision$
 */
public abstract class Aligner
{
    /**
     * Make one.
     * To minimize memory use, seqX ought to be the shorter sequence.
     * Unrecognized characters in the sequence are silently converted to X
     * (which is the ambiguity code we're using to mean "not a base").
     */
    public Aligner( CharSequence seqX, CharSequence seqY, Scorer scorer, boolean lateGapping )
    {
        if ( seqX == null || seqX.length() == 0 )
        {
            throw new IllegalArgumentException("Empty X sequence.");
        }
        if ( seqY == null || seqY.length() == 0 )
        {
            throw new IllegalArgumentException("Empty Y sequence.");
        }
        if ( scorer == null )
        {
            throw new IllegalArgumentException("Scorer is null.");
        }

        mSeqX = seqX;
        mSeqXLen = seqX.length();
        mSeqY = seqY;
        mSeqYLen = seqY.length();
        mScorer = scorer;
        setLateGapping(lateGapping);
    }

    public boolean isLateGapping()
    {
        return mComparator == gLateGapComparator;
    }

    public void setLateGapping( boolean lateGapping )
    {
        mComparator = lateGapping ? gLateGapComparator : gEarlyGapComparator;
    }

    /**
     * The X sequence to align.
     */
    public CharSequence getXSequence()
    {
        return mSeqX;
    }

    /**
     * The Y sequence to align.
     */
    public CharSequence getYSequence()
    {
        return mSeqY;
    }

    /**
     * The scoring matrix used in aligning.
     */
    public Scorer getScorer()
    {
        return mScorer;
    }

    /**
     * This method does the alignment, but does not save the traceback.
     * It is O(NM) in time, and O(N) in space (where N is the length of
     * sequence X and M is the length of sequence Y).
     */
    public abstract Score getScore();

    /**
     * This method does the complete alignment.
     * It is O(NM) in both time and space.  The actual space required (in bytes)
     * is NM/2 + 8N.  (Or half that for self alignments.)
     * So you can't align whole chromosomes, but you can get pretty far.
     */
    public abstract Alignment getAlignment();

    /**
     * Clone this aligner for some other pair of sequences.
     * The resulting aligner will have the same alignment strategy (Global, semi-global, etc.),
     * scoring matrix, late- or early-gapping philosophy, etc.
     */
    public abstract Aligner clone( CharSequence seqX, CharSequence seqY );

    protected Score getScore( Score score, float rowInit, float rowIncr, float colInit, float colIncr )
    {
        // this method implements the recurrences described on p.244 of Gusfield's book.

        int nRows = mSeqYLen;
        int nCols = mSeqXLen;
        int lastCol = nCols - 1;

        float[] maxScores = new float[nCols]; // the max value of the 3 recurrences for a vertical, diagonal, or horizontal movement
        float[] vScores = new float[nCols]; // the value for a vertical move (introduces or extends a gap in X)

        Scorer scorer = mScorer;
        float vGapOpen = scorer.getXGapScorer().getOpenScore(); // cost of moving vertically when starting a gap
        float vGapExtend = scorer.getXGapScorer().getExtendScore(); // cost of moving vertically when extending a gap
        float hGapOpen = scorer.getYGapScorer().getOpenScore(); // cost of moving horizontally when starting a gap
        float hGapExtend = scorer.getYGapScorer().getExtendScore(); // cost of moving horizontally when extending a gap

        int baseX0 = BaseCall.ordinalOf(mSeqX.charAt(0));
        int baseY = BaseCall.ordinalOf(mSeqY.charAt(0));

        // (0,0) cannot extend any gaps
        float dScore = scorer.getScore(baseX0,baseY);
        float hScore = rowInit + hGapOpen;
        float vScore = colInit + vGapOpen;
        float maxScore = Math.max(dScore,Math.max(hScore,vScore));

        vScores[0] = vScore;
        maxScores[0] = maxScore;

        // (1..n,0) cannot extend vGap
        for ( int col = 1; col < nCols; ++col )
        {
            dScore = colInit + scorer.getScore(BaseCall.ordinalOf(mSeqX.charAt(col)),baseY);
            hScore = Math.max(maxScore + hGapOpen,hScore + hGapExtend);
            colInit += colIncr;
            maxScore = Math.max(dScore,Math.max(hScore,colInit + vGapOpen));

            vScores[col] = vScore;
            maxScores[col] = maxScore;
        }

        score.checkRowEnd(maxScore,lastCol,0);

        float maxScorePrevRow;
        for ( int row = 1; row < nRows; ++row )
        {
            baseY = BaseCall.ordinalOf(mSeqY.charAt(row));

            // (0,1..m) can't extend hGap
            dScore = rowInit + scorer.getScore(baseX0,baseY);
            rowInit += rowIncr;
            hScore = rowInit + hGapOpen;
            maxScorePrevRow = maxScores[0];
            vScore = Math.max(maxScorePrevRow + vGapOpen,vScores[0] + vGapExtend);
            maxScore = Math.max(dScore,Math.max(hScore,vScore));

            vScores[0] = vScore;
            maxScores[0] = maxScore;

            // (1..n,1..m) general case
            for ( int col = 1; col < nCols; ++col )
            {
                dScore = maxScorePrevRow + scorer.getScore(BaseCall.ordinalOf(mSeqX.charAt(col)),baseY);
                hScore = Math.max(maxScore + hGapOpen,hScore + hGapExtend);
                maxScorePrevRow = maxScores[col];
                vScore = Math.max(maxScorePrevRow + vGapOpen,vScores[col] + vGapExtend);
                maxScore = Math.max(dScore,Math.max(hScore,vScore));

                vScores[col] = vScore;
                maxScores[col] = maxScore;
            }

            score.checkRowEnd(maxScore,lastCol,row);
        }

        score.checkLastRow(maxScores,nRows-1);

        return score;
    }

    protected LinkedList<Block> getAlignment( Score score, float maxScorePrevColInit, float maxScorePrevColIncr, float maxScorePrevRowInit, float maxScorePrevRowIncr )
    {
        NybbleMatrix traceback = fillTracebackMatrix(score,maxScorePrevColInit,maxScorePrevColIncr,maxScorePrevRowInit,maxScorePrevRowIncr);
        return runTraceback(traceback,score);
    }

    private NybbleMatrix fillTracebackMatrix( Score score, float rowInit, float rowIncr, float colInit, float colIncr )
    {
        int nRows = mSeqYLen;
        int nCols = mSeqXLen;
        int lastCol = nCols - 1;

        NybbleMatrix traceback = new NybbleMatrix(nCols,nRows);
        SerialFiller filler = traceback.new SerialFiller();

        float[] maxScores = new float[nCols]; // the max value of the 3 recurrences for a vertical, diagonal, or horizontal movement
        float[] vScores = new float[nCols]; // the value for a vertical move (introduces or extends a gap in X)

        Scorer scorer = mScorer;
        float vGapOpen = scorer.getXGapScorer().getOpenScore(); // cost of moving vertically when starting a gap
        float vGapExtend = scorer.getXGapScorer().getExtendScore(); // cost of moving vertically when extending a gap
        float hGapOpen = scorer.getYGapScorer().getOpenScore(); // cost of moving horizontally when starting a gap
        float hGapExtend = scorer.getYGapScorer().getExtendScore(); // cost of moving horizontally when extending a gap

        int baseX0 = BaseCall.ordinalOf(mSeqX.charAt(0));
        int baseY = BaseCall.ordinalOf(mSeqY.charAt(0));

        // (0,0) all 3 previous scores are off the edge -- cannot extend any gaps
        float cScore = scorer.getScore(baseX0,baseY);
        float dScore = cScore;
        float hScore = rowInit + hGapOpen;
        float vScore = colInit + vGapOpen;
        int tracebackFlags = TRACEBACK_GAP_OPEN_VERTICAL|TRACEBACK_GAP_OPEN_HORIZONTAL;
        float maxScore = mComparator.test(filler,cScore,tracebackFlags,dScore,hScore,vScore);

        vScores[0] = vScore;
        maxScores[0] = maxScore;

        float openScore;

        // (1..n,0) previous vertical and diagonal scores are off the edge -- cannot extend vGap
        for ( int col = 1; col < nCols; ++col )
        {
            cScore = scorer.getScore(BaseCall.ordinalOf(mSeqX.charAt(col)),baseY);
            dScore = colInit + cScore;
            hScore += hGapExtend;
            tracebackFlags = TRACEBACK_GAP_OPEN_VERTICAL;
            openScore = maxScore + hGapOpen;
            if ( openScore >= hScore )
            {
                hScore = openScore;
                tracebackFlags = TRACEBACK_GAP_OPEN_VERTICAL|TRACEBACK_GAP_OPEN_HORIZONTAL;
            }
            colInit += colIncr;
            vScore = colInit + vGapOpen;
            maxScore = mComparator.test(filler,cScore,tracebackFlags,dScore,hScore,vScore);

            vScores[col] = vScore;
            maxScores[col] = maxScore;
        }

        score.checkRowEnd(maxScore,lastCol,0);

        float maxScorePrevRow;

        for ( int row = 1; row < nRows; ++row )
        {
            baseY = BaseCall.ordinalOf(mSeqY.charAt(row));

            // (0,1..m) previous horizonal score is off the edge -- can't extend hGap
            cScore = scorer.getScore(baseX0,baseY);
            dScore = rowInit + cScore;
            rowInit += rowIncr;
            hScore = rowInit + hGapOpen;
            tracebackFlags = TRACEBACK_GAP_OPEN_HORIZONTAL;
            maxScorePrevRow = maxScores[0];
            vScore = vScores[0] + vGapExtend;
            openScore = maxScorePrevRow + vGapOpen;
            if ( openScore >= vScore )
            {
                vScore = openScore;
                tracebackFlags = TRACEBACK_GAP_OPEN_VERTICAL|TRACEBACK_GAP_OPEN_HORIZONTAL;
            }
            maxScore = mComparator.test(filler,cScore,tracebackFlags,dScore,hScore,vScore);

            vScores[0] = vScore;
            maxScores[0] = maxScore;

            // (1..n,1..m) general case
            for ( int col = 1; col < nCols; ++col )
            {
                cScore = scorer.getScore(BaseCall.ordinalOf(mSeqX.charAt(col)),baseY);
                dScore = maxScorePrevRow + cScore;
                hScore += hGapExtend;
                tracebackFlags = 0;
                openScore = maxScore + hGapOpen;
                if ( openScore >= hScore )
                {
                    hScore = openScore;
                    tracebackFlags = TRACEBACK_GAP_OPEN_HORIZONTAL;
                }
                maxScorePrevRow = maxScores[col];
                vScore = vScores[col] + vGapExtend;
                openScore = maxScorePrevRow + vGapOpen;
                if ( openScore >= vScore )
                {
                    vScore = openScore;
                    tracebackFlags |= TRACEBACK_GAP_OPEN_VERTICAL;
                }
                maxScore = mComparator.test(filler,cScore,tracebackFlags,dScore,hScore,vScore);

                vScores[col] = vScore;
                maxScores[col] = maxScore;
            }

            score.checkRowEnd(maxScore,lastCol,row);
        }

        score.checkLastRow(maxScores,nRows-1);
        filler.done();

        return traceback;
    }

    private LinkedList<Block> runTraceback( NybbleMatrix traceback, Score score )
    {
        LinkedList<Block> blocks = new LinkedList<Block>();
        int xIdx = score.getXBestIndex();
        int yIdx = score.getYBestIndex();
        while ( xIdx >= 0 && yIdx >= 0 )
        {
            int blockEndXIdx = xIdx;
            int blockEndYIdx = yIdx;

            int tracebackVal = traceback.get(xIdx,yIdx);
            int dir = tracebackVal & TRACEBACK_DIRECTION_MASK;
            if ( dir == TRACEBACK_HORIZONTAL )
            {
                while ( (tracebackVal & TRACEBACK_GAP_OPEN_HORIZONTAL) == 0 )
                {
                    tracebackVal = traceback.get(--xIdx,yIdx);
                }

                CharSubSequence xSeq = new CharSubSequence(mSeqX,xIdx--,blockEndXIdx+1);
                CharSubSequence ySeq = new CharSubSequence(mSeqY,yIdx+1,yIdx+1);
                blocks.addFirst(new Block(xSeq,ySeq,0));
            }
            else if ( dir == TRACEBACK_VERTICAL )
            {
                while ( (tracebackVal & TRACEBACK_GAP_OPEN_VERTICAL) == 0 )
                {
                    tracebackVal = traceback.get(xIdx,--yIdx);
                }

                CharSubSequence xSeq = new CharSubSequence(mSeqX,xIdx+1,xIdx+1);
                CharSubSequence ySeq = new CharSubSequence(mSeqY,yIdx--,blockEndYIdx+1);
                blocks.addFirst(new Block(xSeq,ySeq,0));
            }
            else
            {
                int nMatches = 0;
                do
                {
                    if ( dir == TRACEBACK_MATCH )
                    {
                        nMatches += 1;
                    }
                    xIdx -= 1;
                    yIdx -= 1;
                    if ( xIdx < 0 || yIdx < 0 )
                    {
                        break;
                    }
                    dir = traceback.get(xIdx,yIdx) & TRACEBACK_DIRECTION_MASK;
                }
                while ( dir == TRACEBACK_MATCH || dir == TRACEBACK_MISMATCH );

                CharSubSequence xSeq = new CharSubSequence(mSeqX,xIdx+1,blockEndXIdx+1);
                CharSubSequence ySeq = new CharSubSequence(mSeqY,yIdx+1,blockEndYIdx+1);
                blocks.addFirst(new Block(xSeq,ySeq,nMatches));
            }
        }

        if ( xIdx >= 0 )
        {
            CharSubSequence xSeq = new CharSubSequence(mSeqX,0,xIdx+1);
            CharSubSequence ySeq = new CharSubSequence(mSeqY,0,0);
            blocks.addFirst(new Block(xSeq,ySeq,0));
        }
        if ( yIdx >= 0 )
        {
            CharSubSequence xSeq = new CharSubSequence(mSeqX,0,0);
            CharSubSequence ySeq = new CharSubSequence(mSeqY,0,yIdx+1);
            blocks.addFirst(new Block(xSeq,ySeq,0));
        }
        return blocks;
    }

    protected CharSequence mSeqX;
    protected int mSeqXLen;
    protected CharSequence mSeqY;
    protected int mSeqYLen;
    protected Scorer mScorer;
    private ScoreComparator mComparator;

    private static ScoreComparator gEarlyGapComparator = new EarlyGapComparator();
    private static ScoreComparator gLateGapComparator = new LateGapComparator();

    static final int TRACEBACK_MATCH = 0;
    static final int TRACEBACK_HORIZONTAL = 1;
    static final int TRACEBACK_VERTICAL = 2;
    static final int TRACEBACK_MISMATCH = 3;
    static final int TRACEBACK_DIRECTION_MASK = 3;
    static final int TRACEBACK_GAP_OPEN_HORIZONTAL = 4;
    static final int TRACEBACK_GAP_OPEN_VERTICAL = 8;

    interface ScoreComparator
    {
        float test( SerialFiller filler, float cScore, int flags, float dScore, float hScore, float vScore );
    }

    static final class EarlyGapComparator
        implements ScoreComparator
    {
        public final float test( SerialFiller filler, float cScore, int flags, float dScore, float hScore, float vScore )
        {
            float result;

            if ( dScore >= hScore )
            {
                if ( dScore >= vScore )
                {
                    if ( cScore < 0.F )
                    {
                        flags |= TRACEBACK_MISMATCH;
                    }
//                  else
//                  {
//                      val |= TRACEBACK_MATCH;  // don't have to do this because TRACEBACK_MATCH == 0
//                  }
                    result = dScore;
                }
                else if ( hScore >= vScore )
                {
                    flags |= TRACEBACK_HORIZONTAL;
                    result = hScore;
                }
                else
                {
                    flags |= TRACEBACK_VERTICAL;
                    result = vScore;
                }
            }
            else if ( hScore >= vScore )
            {
                flags |= TRACEBACK_HORIZONTAL;
                result = hScore;
            }
            else
            {
                flags |= TRACEBACK_VERTICAL;
                result = vScore;
            }

            filler.set(flags);

            return result;
        }
    }

    static final class LateGapComparator
        implements ScoreComparator
    {
        public final float test( SerialFiller filler, float cScore, int flags, float dScore, float hScore, float vScore )
        {
            float result;

            if ( dScore > hScore )
            {
                if ( dScore > vScore )
                {
                    if ( cScore < 0.F )
                    {
                        flags |= TRACEBACK_MISMATCH;
                    }
//                  else
//                  {
//                      val |= TRACEBACK_MATCH;  // don't have to do this because TRACEBACK_MATCH == 0
//                  }
                    result = dScore;
                }
                else if ( hScore >= vScore )
                {
                    flags |= TRACEBACK_HORIZONTAL;
                    result = hScore;
                }
                else
                {
                    flags |= TRACEBACK_VERTICAL;
                    result = vScore;
                }
            }
            else if ( hScore >= vScore )
            {
                flags |= TRACEBACK_HORIZONTAL;
                result = hScore;
            }
            else
            {
                flags |= TRACEBACK_VERTICAL;
                result = vScore;
            }

            filler.set(flags);

            return result;
        }
    }
}
