/**
 * $Id: SelfAligner.java 73459 2008-10-10 16:10:03Z tsharpe $
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

/**
 * A sequence is aligned against itself.
 * Obviously, the main diagonal aligns with 100% identity, and is uninteresting.
 * Also, the two triangles separated by the main diagonal are mirror images, and
 * so only one of the triangles needs to be computed.  So we just run the local
 * alignment recurrences in the lower-left triangle.  We consider the main diagonal,
 * as well as the left edge, to be "off the edge" of the DP matrix, and give these
 * cells a score of 0.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class SelfAligner
    extends Aligner
{
    public SelfAligner( CharSequence seq )
    {
        super(seq,seq,new Scorer(),false);
    }

    public SelfAligner( CharSequence seq, Scorer scorer )
    {
        super(seq,seq,scorer,false);
    }

    public CharSequence getSequence()
    {
        return mSeqX;
    }

    @Override
    public Score getScore()
    {
        Score score = new Score();
        CharSequence seq = getSequence();
        int nRows = seq.length();
        float[] maxScores = new float[nRows-1]; // the max value of the 3 recurrences for a vertical, diagonal, or horizontal movement
        float[] vScores = new float[nRows-1]; // the value for a vertical move (introduces or extends a gap in X)

        Scorer scorer = getScorer();
        float vGapOpen = scorer.getXGapScorer().getOpenScore(); // cost of moving vertically when starting a gap
        float vGapExtend = scorer.getXGapScorer().getExtendScore(); // cost of moving vertically when extending a gap
        float hGapOpen = scorer.getYGapScorer().getOpenScore(); // cost of moving horizontally when starting a gap
        float hGapExtend = scorer.getYGapScorer().getExtendScore(); // cost of moving horizontally when extending a gap

        int baseX0 = BaseCall.ordinalOf(seq.charAt(0));
        int baseY = BaseCall.ordinalOf(seq.charAt(1));

        // (0,1) is off the edge in all 3 directions -- can't extend any gaps
        float maxScore = Math.max(scorer.getScore(baseX0,baseY),Math.max(hGapOpen,vGapOpen));
        vScores[0] = vGapOpen;
        maxScores[0] = maxScore;
        score.checkScore(maxScore, 0, 1);

        for ( int row = 2; row < nRows; ++row )
        {
            baseY = BaseCall.ordinalOf(seq.charAt(row));

            // (0,2..nRows) is on the left edge so previous column and diagonal scores are off the edge -- can't extend horizontal gap
            float dScore = scorer.getScore(baseX0,baseY);
            float hScore = hGapOpen;
            float maxScorePrevRow = maxScores[0];
            float vScore = Math.max(vScores[0] + vGapExtend,maxScorePrevRow + vGapOpen);
            maxScore = Math.max(dScore,Math.max(hScore,vScore));
            vScores[0] = vScore;
            maxScores[0] = maxScore;
            score.checkScore(maxScore, 0, row);

            // (1..row-1,2..nRows) is the general case
            int maxCol = row - 1;
            for ( int col = 1; col < maxCol; ++col )
            {
                dScore = maxScorePrevRow + scorer.getScore(BaseCall.ordinalOf(seq.charAt(col)),baseY);
                hScore = Math.max(hScore + hGapExtend,maxScore + hGapOpen);
                maxScorePrevRow = maxScores[col];
                vScore = Math.max(vScores[col] + vGapExtend,maxScorePrevRow + vGapOpen);
                maxScore = Math.max(dScore,Math.max(hScore,vScore));
                vScores[col] = vScore;
                maxScores[col] = maxScore;
                score.checkScore(maxScore, col, row);
            }

            // (row-1,2..nRows) is below the main diagonal -- can't extend vertical gap
            dScore = maxScorePrevRow + scorer.getScore(BaseCall.ordinalOf(seq.charAt(maxCol)),baseY);
            hScore = Math.max(hScore + hGapExtend,maxScore + hGapOpen);
            maxScore = Math.max(dScore,Math.max(hScore,vGapOpen));
            vScores[maxCol] = vGapOpen;
            maxScores[maxCol] = maxScore;
            score.checkScore(maxScore, maxCol, row);
        }
        return score;
    }

    @Override
    public Alignment getAlignment()
    {
        Score score = new Score();
        NybbleDiagMatrix traceback = calcTraceback(score);
        return runTraceback(traceback,score);
    }

    @Override
    public Aligner clone( CharSequence seqX, CharSequence seqY )
    {
        if ( seqX != seqY )
        {
            throw new IllegalArgumentException("seqX and seqY must be identical for a self-alignment");
        }
        return new SelfAligner(seqX,getScorer());
    }

    private NybbleDiagMatrix calcTraceback( Score score )
    {
        CharSequence seq = getSequence();
        int nRows = seq.length();
        NybbleDiagMatrix traceback = new NybbleDiagMatrix(nRows);
        float[] maxScores = new float[nRows-1]; // the max value of the 3 recurrences for a vertical, diagonal, or horizontal movement
        float[] vScores = new float[nRows-1]; // the value for a vertical move (introduces or extends a gap in X)

        Scorer scorer = getScorer();
        float vGapOpen = scorer.getXGapScorer().getOpenScore(); // cost of moving vertically when starting a gap
        float vGapExtend = scorer.getXGapScorer().getExtendScore(); // cost of moving vertically when extending a gap

        float hGapOpen = scorer.getYGapScorer().getOpenScore(); // cost of moving horizontally when starting a gap
        float hGapExtend = scorer.getYGapScorer().getExtendScore(); // cost of moving horizontally when extending a gap

        int baseX0 = BaseCall.ordinalOf(seq.charAt(0));
        int baseY = BaseCall.ordinalOf(seq.charAt(1));

        // (0,1) is off the edge in all 3 directions -- can't extend any gaps
        int tracebackVal = TRACEBACK_GAP_OPEN_HORIZONTAL|TRACEBACK_GAP_OPEN_VERTICAL;
        float maxScore = setTraceback(traceback,0,1,tracebackVal,scorer.getScore(baseX0,baseY),hGapOpen,vGapOpen);
        vScores[0] = vGapOpen;
        maxScores[0] = maxScore;
        score.checkScore(maxScore, 0, 1);

        for ( int row = 2; row < nRows; ++row )
        {
            baseY = BaseCall.ordinalOf(seq.charAt(row));

            // (0,2..nRows) is on the left edge so previous column and diagonal scores are off the edge -- can't extend horizontal gap
            float dScore = scorer.getScore(baseX0,baseY);
            float hScore = hGapOpen;
            tracebackVal = TRACEBACK_GAP_OPEN_HORIZONTAL;
            float vScore = vScores[0] + vGapExtend;
            float maxScorePrevRow = maxScores[0];
            float openScore = maxScorePrevRow + vGapOpen;
            if ( openScore >= vScore )
            {
                vScore = openScore;
                tracebackVal = TRACEBACK_GAP_OPEN_VERTICAL | TRACEBACK_GAP_OPEN_HORIZONTAL;
            }
            maxScore = setTraceback(traceback,0,row,tracebackVal,dScore,hScore,vScore);
            vScores[0] = vScore;
            maxScores[0] = maxScore;
            score.checkScore(maxScore, 0, row);

            // (1..row-1,2..nRows) is the general case
            int maxCol = row - 1;
            for ( int col = 1; col < maxCol; ++col )
            {
                dScore = maxScorePrevRow + scorer.getScore(BaseCall.ordinalOf(seq.charAt(col)),baseY);
                hScore += hGapExtend;
                openScore = maxScore + hGapOpen;
                tracebackVal = 0;
                if ( openScore >= hScore )
                {
                    hScore = openScore;
                    tracebackVal = TRACEBACK_GAP_OPEN_HORIZONTAL;
                }
                vScore = vScores[col] + vGapExtend;
                maxScorePrevRow = maxScores[col];
                openScore = maxScorePrevRow + vGapOpen;
                if ( openScore >= vScore )
                {
                    vScore = openScore;
                    tracebackVal |= TRACEBACK_GAP_OPEN_VERTICAL;
                }
                maxScore = setTraceback(traceback,col,row,tracebackVal,dScore,hScore,vScore);
                vScores[col] = vScore;
                maxScores[col] = maxScore;
                score.checkScore(maxScore, col, row);
            }

            // (row-1,2..nRows) is below the main diagonal -- can't extend vertical gap
            dScore = maxScorePrevRow + scorer.getScore(BaseCall.ordinalOf(seq.charAt(maxCol)),baseY);
            tracebackVal = TRACEBACK_GAP_OPEN_VERTICAL;
            hScore += hGapExtend;
            openScore = maxScore + hGapOpen;
            if ( openScore >= hScore )
            {
                hScore = openScore;
                tracebackVal |= TRACEBACK_GAP_OPEN_HORIZONTAL;
            }
            maxScore = setTraceback(traceback,maxCol,row,tracebackVal,dScore,hScore,vGapOpen);
            vScores[maxCol] = vGapOpen;
            maxScores[maxCol] = maxScore;
            score.checkScore(maxScore, maxCol, row);
        }

        return traceback;
    }

    private Alignment runTraceback( NybbleDiagMatrix traceback, Score score )
    {
        LinkedList<Block> blocks = new LinkedList<Block>();
        Scorer scorer = getScorer();
        int xIdx = score.getXBestIndex();
        int yIdx = score.getYBestIndex();
        while ( xIdx >= 0 && yIdx >= 0 && xIdx < yIdx )
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

                CharSubSequence seqX = new CharSubSequence(mSeqX,xIdx--,blockEndXIdx+1);
                CharSubSequence seqY = new CharSubSequence(mSeqY,yIdx+1,yIdx+1);
                blocks.addFirst(new Block(seqX,seqY,0));
            }
            else if ( dir == TRACEBACK_VERTICAL )
            {
                while ( (tracebackVal & TRACEBACK_GAP_OPEN_VERTICAL) == 0 )
                {
                    tracebackVal = traceback.get(xIdx,--yIdx);
                }

                CharSubSequence seqX = new CharSubSequence(mSeqX,xIdx+1,xIdx+1);
                CharSubSequence seqY = new CharSubSequence(mSeqY,yIdx--,blockEndYIdx+1);
                blocks.addFirst(new Block(seqX,seqY,0));
            }
            else if ( dir == TRACEBACK_DIAGONAL )
            {
                int nMatches = 0;
                do
                {
                    int xVal = BaseCall.ordinalOf(getSequence().charAt(xIdx));
                    int yVal = BaseCall.ordinalOf(getYSequence().charAt(yIdx));
                    if ( scorer.getScore(xVal,yVal) >= 0.F )
                    {
                        nMatches += 1;
                    }
                    xIdx -= 1;
                    yIdx -= 1;
                }
                while ( xIdx >= 0 && yIdx >= 0 && (traceback.get(xIdx,yIdx) & TRACEBACK_DIRECTION_MASK) == TRACEBACK_DIAGONAL );

                CharSubSequence seqX = new CharSubSequence(mSeqX,xIdx+1,blockEndXIdx+1);
                CharSubSequence seqY = new CharSubSequence(mSeqY,yIdx+1,blockEndYIdx+1);
                blocks.addFirst(new Block(seqX,seqY,nMatches));
            }
            else
            {
                assert( dir == TRACEBACK_STOP );
                break;
            }
        }

        CharSubSequence seqX = new CharSubSequence(mSeqX,0,mSeqX.length());
        CharSubSequence seqY = new CharSubSequence(mSeqY,0,mSeqY.length());
        return new Alignment(score.getScore(),blocks,seqX,seqY);
    }

    private static float setTraceback( NybbleDiagMatrix traceback, int col, int row, int val, float dScore, float hScore, float vScore )
    {
        float result;
        int dir;
        if ( dScore >= hScore )
        {
            if ( dScore >= vScore )
            {
                dir = TRACEBACK_DIAGONAL;
                result = dScore;
            }
            else if ( hScore >= vScore )
            {
                dir = TRACEBACK_HORIZONTAL;
                result = hScore;
            }
            else
            {
                dir = TRACEBACK_VERTICAL;
                result = vScore;
            }
        }
        else if ( hScore >= vScore )
        {
            dir = TRACEBACK_HORIZONTAL;
            result = hScore;
        }
        else
        {
            dir = TRACEBACK_VERTICAL;
            result = vScore;
        }

        if ( result < 0.F )
        {
            dir = TRACEBACK_STOP;
            result = 0.F;
        }

        traceback.set(col,row,val|dir);

        return result;
    }

    private static final int TRACEBACK_STOP = 0;
    private static final int TRACEBACK_DIAGONAL = 3;
}
