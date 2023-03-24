/*
 * $Id: SuffixTree.java 73459 2008-10-10 16:10:03Z tsharpe $
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
import java.util.List;

/**
 * Implementation of a suffix tree for the alphabet {A,C,G,T}.
 * This is a data structure for exact (sub)string matching.
 *
 * @author tsharpe
 * @version $Revision: 48897 $
 */
public class SuffixTree
{
    /**
     * Make one.
     */
    public SuffixTree()
    {
        mRoot = new SuffixTreeInternalNode(null,0,0);
    }

    /**
     * Make one from a single sequence.
     * You can always add more later.
     * This is a convenience method equivalent to making a SuffixTree, then adding a sequence.
     */
    public SuffixTree( CharSequence str )
    {
        this();
        add(str);
    }

    /**
     * Add a string to the tree.
     * The sequence is broken into runs of [ACGT]+, and each of the runs is added to the tree.
     * @param str The string to add.
     */
    public void add( CharSequence str )
    {
        if ( str == null )
        {
            throw new IllegalArgumentException("Sequence to be added to suffix tree cannot be null.");
        }

        int nnn = str.length();
        for ( int idx = 0; idx < nnn; ++idx )
        {
            // find a valid character
            if ( "ACGT".indexOf(Character.toUpperCase(str.charAt(idx))) != -1 )
            {
                int idx2 = idx;
                while ( ++idx2 < nnn )
                {
                    // break on an invalid character
                    if ( "ACGT".indexOf(Character.toUpperCase(str.charAt(idx2))) == -1 )
                        break;
                }

                CharSequence goodBit = str;    // assume we're adding the whole string
                if ( idx != 0 || idx2 != nnn ) // if we're not adding the whole string
                {
                    goodBit = new CharSubSequence(str,idx,idx2); // take the right piece
                }
                buildImplicitTree(goodBit);
                markTree(goodBit);
                idx = idx2; // outer loop skips ahead to next unexamined character
            }
        }
    }

    /**
     * Returns a list of SuffixTreeMatches.  The matches represent the set of target strings
     * having the longest region of identity with the probe string.  (It's a list because there
     * may be more than one string in the database with an identity of the maximum size.)
     * @param probe The probe sequence.
     * @param minMatchLength The least lengthy match you'll find interesting.
     * @return A list of SuffixTreeMatches.
     */
    public List<SuffixTreeMatch> bestMatches( CharSequence probe, int minMatchLength )
    {
        if ( probe == null )
            throw new IllegalArgumentException("Probe cannot be null.");

        SuffixTreeBestMatchListener listener = new SuffixTreeBestMatchListener(minMatchLength);

        match( probe, listener );

        return listener.getBestMatches();
    }

    /**
     * Find strings in the database that match the entire probe.  The return is in the form of
     * a list of CharSubSequences that contain the probe sequence in its entirety.
     * @param probe The probe sequence.
     * @return A list of CharSubSequences.
     */
    public List<CharSubSequence> probeMatches( CharSequence probe )
    {
        if ( probe == null )
            throw new IllegalArgumentException("Probe cannot be null.");

        SuffixTreeBestMatchListener listener = new SuffixTreeBestMatchListener(probe.length());

        match( probe, listener );

        return listener.getTargets();
    }

    /**
     * Here's how you search the tree.  Provide a probe string and an object to
     * receive a stream of match events.  Your object can save the best match,
     * sort all matches by length, ignore matches shorter than some minimum length,
     * or whatever.
     *
     * @param probe The string to look for.
     * @param listener An object that will receive match events.
     */
    public void match( CharSequence probe, SuffixTreeMatchListener listener )
    {
        if ( probe == null )
            throw new IllegalArgumentException("Probe cannot be null.");
        if ( listener == null )
            throw new IllegalArgumentException("Listener cannot be null.");

        SuffixTreePosition pos = new SuffixTreePosition(mRoot);
        int strLen = probe.length();
        int jjj = 0;
        for ( int iii = 0; iii < strLen; ++iii )
        {
            if ( jjj < iii )
                jjj = iii;

            while ( jjj < strLen && pos.step(probe.charAt(jjj)) )
                jjj += 1;

            if ( !listener.match(probe,iii,jjj,pos) )
                break;

            pos.toSuffix();
            pos.skip(probe,iii+1,jjj);
        }
    }

    /**
     * Print the tree.
     */
    public void printTree()
    {
        for ( char child : CHILDREN )
        {
            printNode(mRoot.getChild(child),"");
        }
    }

    /**
     * Check the structure of the tree.  (For debugging.)
     */
    public void validate()
    {
        for ( char child : CHILDREN )
        {
            validate(mRoot.getChild(child),"");
        }
    }

    /**
     * Returns the length of the longest prefix of a probe that occurs anywhere within any of the strings in this tree.
     * E.g., if the tree contains "ACGTACGT" and "TGTATGCA" the probe "TACA" would have a prefix match length of 3,
     * since the prefix "TAC" occurs (at the fourth character of the 1st string), but "TACA" doesn't appear anywhere
     * in the tree.
     *
     * @param str The probe sequence.
     * @return The length of the longest prefix of the probe that can be found in the tree.
     */
    public int prefixMatchLength( CharSequence str )
    {
        return new SuffixTreePosition(mRoot).walk(str);
    }

    /*
     * This implements Ukkonen's algorithm for building an implicit suffix tree in linear time.
     * See <i>Algorithms on Strings, Trees, and Sequences</i>
     * by Dan Gusfield, Cambridge Univ. Press, reprinted in 1999, pp. 94-107.
     */
    private void buildImplicitTree( CharSequence str )
    {
        char chr = 'X';
        char chr2;
        SuffixTreeNode newInternalNode = null;
        int strLen = str.length();
        int minExIdx = 0;

        for ( int phaseIdx = prefixMatchLength(str); phaseIdx < strLen; ++phaseIdx ) // phase loop
        {
            int strPos = minExIdx;
            int nToMatch = phaseIdx - minExIdx;

            SuffixTreeNode prevNode = SuffixTreeDummyNode.SINGLETON;
            SuffixTreeNode node = mRoot;
            int labLen = 0;

            for ( int exIdx = minExIdx; exIdx <= phaseIdx; ++exIdx ) // extension loop
            {
                // skip/count to find str[exIdx..phaseIdx] (trick 1)
                while ( nToMatch > labLen )
                {
                    strPos += labLen;
                    nToMatch -= labLen;
                    prevNode = node;
                    node = node.getChild((chr = str.charAt(strPos)));
                    labLen = node.getLabelLength();
                    if ( node.isLeaf() && node.getString() == str )
                    {   // this is our version of trick 3 -- we build a leaf at full length, but adjust
                        // its effective length to keep it bounded for the phase we're executing.
                        labLen = phaseIdx - node.getLabelStart();
                    }
                }

                // extend to include str[phaseIdx]
                if ( nToMatch != labLen ) // we're midway through an edge
                {
                    chr2 = Character.toUpperCase(str.charAt(phaseIdx));
                    if ( chr2 != Character.toUpperCase(node.charAt(nToMatch)) ) // Rule 2a -- ended partway through an edge
                    {                                                           // and the next character doesn't match
                        node = node.split(nToMatch);
                        prevNode.setChild(chr,node);
                        node.setChild(chr2,new SuffixTreeLeafNode(str,phaseIdx));
                        if ( newInternalNode != null )
                        {
                            newInternalNode.setSuffixLink(node);
                        }
                        newInternalNode = node;
                    }
                    else // Rule 3 -- found it in the tree already (within some edge)
                    {
                        minExIdx = exIdx;
                        break; // (trick 2)
                    }
                }
                else if ( node.isLeaf() ) // Rule 1 -- ended at a leaf-end
                {
                    if ( node.getString() != str ) // this is a little something extra, since we're building generalized trees
                    {                              // we may have to swizzle the current leaf if it's on some other string
                        node = node.internalize(nToMatch,str,strPos+labLen);
                        prevNode.setChild(chr,node);
                        if ( newInternalNode != null )
                        {
                            newInternalNode.setSuffixLink(node);
                        }
                        newInternalNode = node;
                    }
                }
                else
                {
                    if ( newInternalNode != null )
                    {
                        newInternalNode.setSuffixLink(node);
                        newInternalNode = null;
                    }

                    if ( node.getChild((chr2 = str.charAt(phaseIdx))) == null ) // Rule 2b -- at end of edge, but there's no
                    {                                                           // path that continues from there
                        node.setChild(chr2,new SuffixTreeLeafNode(str,phaseIdx));
                    }
                    else // Rule 3 -- found it in the tree already (at the end of some edge)
                    {
                        minExIdx = exIdx;
                        break; // (trick 2)
                    }
                }

                if ( exIdx == phaseIdx )
                {
                    break;
                }

                // use suffix links to find starting point for next extension
                if ( (node = node.getSuffixLink()) != null )
                {
                    strPos += nToMatch;
                }
                else if ( prevNode.isRoot() )
                {
                    node = prevNode;
                    strPos = exIdx + 1;
                }
                else
                {
                    node = prevNode.getSuffixLink();
                }
                labLen = node.getLabelLength();
                strPos -= labLen;
                nToMatch = phaseIdx - strPos;
            }
        }
    }

    /*
     * This is like an extra, final phase, extending the tree to encompass a fictitious terminator character.
     * (Gusfield, op. cit., p. 107.)
     * (I.e., it transforms an implicit suffix tree for the string into an explicit one.)
     * While we're doing this, we mark each node that terminates a suffix of this string as such.
     * Note that we're making trees sort of halfway between implicit and explicit trees:  Each node has a list
     * of strings that terminate there (including interior nodes), and we don't actually put the terminator
     * character into the leaves.
     */
    private void markTree( CharSequence str )
    {
        int strLen = str.length();
        int strPos = 0;
        int nToMatch = strLen;

        SuffixTreeNode prevNode = SuffixTreeDummyNode.SINGLETON;
        SuffixTreeNode node = mRoot;
        int labLen = 0;

        char chr = 'X';
        SuffixTreeNode newInternalNode = null;

        for ( int exIdx = 0; exIdx < strLen; ++exIdx )
        {
            while ( nToMatch > labLen )
            {
                strPos += labLen;
                nToMatch -= labLen;
                prevNode = node;
                node = node.getChild((chr = str.charAt(strPos)));
                labLen = node.getLabelLength();
            }

            if ( nToMatch != labLen )
            {
                node = node.split(nToMatch);
                prevNode.setChild(chr,node);
                if ( newInternalNode != null )
                {
                    newInternalNode.setSuffixLink(node);
                }
                newInternalNode = node;
            }
            else if ( newInternalNode != null )
            {
                newInternalNode.setSuffixLink(node);
                newInternalNode = null;
            }

            node.addString(str);

            node = node.getSuffixLink();
            if ( node != null )
            {
                strPos += nToMatch;
            }
            else if ( prevNode.isRoot() )
            {
                node = prevNode;
                strPos = exIdx + 1;
            }
            else
            {
                node = prevNode.getSuffixLink();
            }
            labLen = node.getLabelLength();
            strPos -= labLen;
            nToMatch = strLen - strPos;
        }

        if ( newInternalNode != null )
        {
            newInternalNode.setSuffixLink(mRoot);
        }
    }

    private void validate( SuffixTreeNode node, String prefix )
    {
        if ( node != null )
        {
            prefix += node.getLabel();

            if ( !node.isLeaf() )
            {
                SuffixTreeNode suffixNode = node.getSuffixLink();
                if ( suffixNode == null )
                {
                    System.err.println("Internal node " + prefix + " has no suffix link.");
                }
                else
                {
                    CharSequence suffix = prefix.substring(1);
                    int strPos = 0;
                    SuffixTreeNode searchNode = mRoot;
                    int labLen = 0;
                    int nToMatch = suffix.length();
                    while ( nToMatch > labLen )
                    {
                        strPos += labLen;
                        nToMatch -= labLen;
                        searchNode = searchNode.getChild(suffix.charAt(strPos));
                        labLen = searchNode.getLabelLength();
                    }
                    if ( nToMatch != labLen )
                    {
                        System.err.println("Suffix link has wrong shape.");
                    }
                    else if ( searchNode != suffixNode )
                    {
                        System.err.println("Wrong suffix link.");
                    }

                }
            }

            Iterator<CharSequence> itr = node.getStringIterator();
            if ( itr.hasNext() )
            {
                int prefixLen = prefix.length();
                do
                {
                    CharSequence str = itr.next();
                    int len = str.length();
                    CharSequence match = str.subSequence(len-prefixLen,len);
                    if ( !prefix.equals(match) )
                        System.err.print(prefix + " != " + match);
                }
                while ( itr.hasNext() );
            }

            for ( char child : CHILDREN )
            {
                validate(node.getChild(child),prefix);
            }
        }
    }

    private void printNode( SuffixTreeNode node, String prefix )
    {
        if ( node != null )
        {
            prefix += node.getLabel();

            Iterator<CharSequence> itr = node.getStringIterator();
            if ( itr.hasNext() )
            {
                System.out.print(prefix);
                System.out.print(" is a suffix of");
                do
                {
                    CharSequence str = itr.next();
                    System.out.print(' ');
                    System.out.print( str );
                }
                while ( itr.hasNext() );
                System.out.println();
            }

            for ( char child : CHILDREN )
            {
                printNode(node.getChild(child),prefix);
            }
        }
    }

    protected SuffixTreeNode mRoot;
    protected static final char[] CHILDREN = {'A','C','G','T'};
}
