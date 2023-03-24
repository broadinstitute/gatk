/*
 * $Id: SuffixTreePosition.java 73459 2008-10-10 16:10:03Z tsharpe $
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

import java.util.Objects;

/**
 * Tracks progress while walking a suffix tree.
 *
 * @author tsharpe
 * @version $Revision: 48897 $
 */
public class SuffixTreePosition
    implements Cloneable
{
    public SuffixTreePosition( SuffixTreeNode root )
    {
        mRoot = root;
        init();
    }

    public int walk( CharSequence str )
    {
        init();
        int strLen = str.length();
        for ( int iii = 0; iii < strLen; ++iii )
            if ( !step(str.charAt(iii)) )
                break;

        return mPathLen + mIdx;
    }

    public boolean step( char chr )
    {
        boolean result = false;

        int labLen = mNode.getLabelLength();
        if ( labLen == mIdx )
        {
            SuffixTreeNode node = mNode.getChild(chr);
            if ( node != null )
            {
                mPrevNode = mNode;
                mNode = node;
                mPathLen += labLen;
                mIdx = 1;
                result = true;
            }
        }
        else if ( Character.toUpperCase(mNode.charAt(mIdx)) == Character.toUpperCase(chr) )
        {
            mIdx += 1;
            result = true;
        }

        return result;
    }

    public void skip( CharSequence str, int start, int stop )
    {
        int depth = stop - start;
        if ( depth < 0 )
        {
            init();
        }
        else
        {
            int labLen;
            while (  mPathLen + (labLen = mNode.getLabelLength()) < depth )
            {
                mPathLen += labLen;
                mPrevNode = mNode;
                mNode = mNode.getChild(str.charAt(start+mPathLen));
            }
            mIdx = depth - mPathLen;
        }
    }

    public void toSuffix()
    {
        SuffixTreeNode node = mPrevNode.getSuffixLink();
        if ( node != null && node != mRoot )
        {
            mPrevNode = mRoot;
            mNode = node;
            mPathLen -= node.getLabelLength() + 1;
            mIdx = 0;
        }
        else
        {
            init();
        }
    }

    public int getMatchLen()
    {
        return mPathLen + mIdx;
    }

    public int getUnmatchedSuffixLen()
    {
        return mNode.getLabelLength() - mIdx;
    }

    public SuffixTreeNode getNode()
    {
        return mNode;
    }

    @Override
    public boolean equals( Object obj )
    {
        boolean result = false;

        if ( this == obj )
        {
            result = true;
        }
        else if ( obj instanceof SuffixTreePosition )
        {
            SuffixTreePosition that = (SuffixTreePosition)obj;
            if ( this.mRoot == that.mRoot &&
                 this.mNode == that.mNode &&
                 this.mPathLen == that.mPathLen &&
                 this.mIdx == that.mIdx )
            {
                result = true;
            }
        }

        return result;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mRoot, mNode, mPathLen, mIdx);
    }

    @Override
    public Object clone()
    {
        try
        {
            return super.clone();
        }
        catch ( CloneNotSupportedException ignore )
        {
            throw new IllegalStateException("Can't clone SuffixTreePosition.");
        }
    }

    private void init()
    {
        mPrevNode = mRoot;
        mNode = mRoot;
        mPathLen = 0;
        mIdx = 0;
    }

    /*
    private int findDepth( SuffixTreeNode target, SuffixTreeNode node, int depth )
    {
        if ( node == target )
            return depth;

        int result = -1;
        if ( node != null )
        {
            depth += node.getLabelLength();
            result = findDepth( target, node.getChild('A'), depth );
            if ( result == -1 )
                result = findDepth( target, node.getChild('C'), depth );
            if ( result == -1 )
                result = findDepth( target, node.getChild('G'), depth );
            if ( result == -1 )
                result = findDepth( target, node.getChild('T'), depth );
        }

        return result;
    }
    */

    private SuffixTreeNode mRoot;
    private SuffixTreeNode mNode;
    private SuffixTreeNode mPrevNode;
    private int mPathLen;
    private int mIdx;
}
