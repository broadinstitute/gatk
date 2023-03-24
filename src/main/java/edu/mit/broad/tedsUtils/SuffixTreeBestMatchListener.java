/*
 * $Id: SuffixTreeBestMatchListener.java 73459 2008-10-10 16:10:03Z tsharpe $
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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * A match listener that finds the longest match.
 *
 * @author tsharpe
 * @version $Revision: 28397 $
 */
public class SuffixTreeBestMatchListener
    implements SuffixTreeMatchListener
{
    public SuffixTreeBestMatchListener( int minMatchLen )
    {
        mMatchLen = minMatchLen;
    }

    /**
     * listener method
     */
    public boolean match( CharSequence probe, int start, int end, SuffixTreePosition pos )
    {
        int matchLen = pos.getMatchLen();

        if ( matchLen > mMatchLen  )
        {
            if ( mMatches != null )
                mMatches.clear();
            mMatchLen = pos.getMatchLen();
        }

        if ( matchLen == mMatchLen )
        {
            if ( mMatches == null )
                mMatches = new ArrayList<Object>();
            mMatches.add( new CharSubSequence(probe,start,end) );
            mMatches.add( pos.clone() );
        }

        return true;
    }

    public boolean hasMatch()
    {
        return mMatches != null;
    }

    public List<SuffixTreeMatch> getBestMatches()
    {
        List<SuffixTreeMatch> list = new ArrayList<SuffixTreeMatch>();

        if ( mMatches != null )
        {
            Iterator<Object> itr = mMatches.iterator();
            while ( itr.hasNext() )
            {
                CharSubSequence probe = (CharSubSequence)itr.next();
                SuffixTreePosition pos = (SuffixTreePosition)itr.next();
                buildMatchList( list, probe, pos.getNode(), pos.getUnmatchedSuffixLen() );
            }
        }

        return list;
    }

    public List<CharSubSequence> getTargets()
    {
        List<CharSubSequence> list = new ArrayList<CharSubSequence>();

        if ( mMatches != null )
        {
            Iterator<Object> itr = mMatches.iterator();
            while ( itr.hasNext() )
            {
                CharSubSequence probe = (CharSubSequence)itr.next();
                SuffixTreePosition pos = (SuffixTreePosition)itr.next();
                buildTargetList( list, probe.length(), pos.getNode(), pos.getUnmatchedSuffixLen() );
            }
        }

        return list;
    }

    private static void buildMatchList( List<SuffixTreeMatch> list, CharSubSequence probe, SuffixTreeNode node, int unmatchedSuffixLen )
    {
        Iterator<CharSequence> itr = node.getStringIterator();
        while ( itr.hasNext() )
        {
            CharSequence target = itr.next();
            int end = target.length() - unmatchedSuffixLen;
            int start = end - probe.length();
            list.add( new SuffixTreeMatch(probe,new CharSubSequence(target,start,end)) );
        }

        SuffixTreeNode child;
        if ( (child = node.getChild('A')) != null )
            buildMatchList( list, probe, child, unmatchedSuffixLen+child.getLabelLength() );
        if ( (child = node.getChild('C')) != null )
            buildMatchList( list, probe, child, unmatchedSuffixLen+child.getLabelLength() );
        if ( (child = node.getChild('G')) != null )
            buildMatchList( list, probe, child, unmatchedSuffixLen+child.getLabelLength() );
        if ( (child = node.getChild('T')) != null )
            buildMatchList( list, probe, child, unmatchedSuffixLen+child.getLabelLength() );
    }

    private static void buildTargetList( List<CharSubSequence> list, int probeLen, SuffixTreeNode node, int unmatchedSuffixLen )
    {
        Iterator<CharSequence> itr = node.getStringIterator();
        while ( itr.hasNext() )
        {
            CharSequence target = itr.next();
            int end = target.length() - unmatchedSuffixLen;
            int start = end - probeLen;
            list.add( new CharSubSequence(target,start,end) );
        }

        SuffixTreeNode child;
        if ( (child = node.getChild('A')) != null )
            buildTargetList( list, probeLen, child, unmatchedSuffixLen+child.getLabelLength() );
        if ( (child = node.getChild('C')) != null )
            buildTargetList( list, probeLen, child, unmatchedSuffixLen+child.getLabelLength() );
        if ( (child = node.getChild('G')) != null )
            buildTargetList( list, probeLen, child, unmatchedSuffixLen+child.getLabelLength() );
        if ( (child = node.getChild('T')) != null )
            buildTargetList( list, probeLen, child, unmatchedSuffixLen+child.getLabelLength() );
    }

    private List<Object> mMatches;
    private int mMatchLen;
}
