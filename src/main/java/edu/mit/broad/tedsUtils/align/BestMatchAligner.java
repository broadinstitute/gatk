/*
 * $Id: BestMatchAligner.java 74578 2008-10-30 21:24:57Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2008 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils.align;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import edu.mit.broad.tedsUtils.CharSubSequence;
import edu.mit.broad.tedsUtils.SuffixTree;
import edu.mit.broad.tedsUtils.SuffixTreeNode;
import edu.mit.broad.tedsUtils.SuffixTreePosition;

/**
 * Finds best matching sequence, and aligns to it.<br>
 * This class is useful, for example, when you have a pile of reads, and a modest
 * number (suffix trees can get big quickly) of amplicons, and need to identify
 * which amplicon a read belongs to and align it against that amplicon.
 *
 * <p>You give it a set of target sequences, and an Aligner to use as an exemplar.
 * All target sequences having an exact match of a specified length are aligned
 * to the probe sequence, and the alignment having the greatest number of matching
 * bases is returned.
 *
 * <p>If the directionality of the probes is not known, you may wish to supply a set
 * of targets that includes both a reference sequence and its reverse complement.
 * If the directionality is known, it would be more efficient just to reverse
 * complement the reverse-direction probes.
 *
 * @author tsharpe
 */
public class BestMatchAligner
{
    /**
     * Make one.
     * @param targets The sequences you want to align to.
     * @param revCompTargets The reverse-complement of those sequences, or null.
     * @param exemplar How you want to do the alignment.
     * @param minMatchLen The minimum exact-match length that will trigger a trial alignment.
     */
    public BestMatchAligner( Collection<? extends CharSequence> targets, Collection<? extends CharSequence> revCompTargets, Aligner exemplar, int minMatchLen )
    {
        if ( targets == null )
        {
            throw new IllegalArgumentException("set of targets may not be null");
        }
        if ( exemplar == null )
        {
            throw new IllegalArgumentException("aligner exemplar may not be null");
        }
        if ( minMatchLen < 7 )
        {
            throw new IllegalArgumentException("minMatchLength too small -- you might as well just align to each target");
        }
        mRevCompTargets = revCompTargets;
        mSuffixTree = new XTree(minMatchLen);
        for ( CharSequence target : targets )
        {
            mSuffixTree.add(target);
        }
        if ( revCompTargets != null )
        {
            for ( CharSequence target : revCompTargets )
            {
                mSuffixTree.add(target);
            }
        }
        mSuffixTree.decorate();
        mExemplar = exemplar;
    }

    /**
     * Find the best alignment.  May return null.
     * The alignment of reverse-complemented sequences will be done with late-gapping
     * set to the opposite sense of the late-gapping flag set in the exemplar.
     */
    public Alignment getAlignment( CharSequence probe )
    {
        if ( probe == null )
        {
            throw new IllegalArgumentException("probe cannot be null");
        }

        Alignment bestAlignment = null;
        int bestCount = 0;

        Set<CharSequence> matches = mSuffixTree.match(probe);
        for ( CharSequence seq : matches )
        {
            Aligner aligner = mExemplar.clone(probe,seq);
            if ( mRevCompTargets.contains(seq) )
            {
                aligner.setLateGapping(!aligner.isLateGapping());
            }
            Alignment alignment = aligner.getAlignment();
            if ( alignment.getNMatches() > bestCount )
            {
                bestAlignment = alignment;
                bestCount = alignment.getNMatches();
            }
        }

        return bestAlignment;
    }

    private Collection<? extends CharSequence> mRevCompTargets;
    private XTree mSuffixTree;
    private Aligner mExemplar;

    /* An extended suffix tree.
     * It decorates each node with the set of target strings that might
     * be matched at a specified minimum length.
     * Each node with a path length greater than the minimum has the set
     * of all strings matched by this node or any of its children.  All this
     * pre-computation is still O(n), I believe, and lets you get the complete
     * set of possible targets for any probe in O(m) time. 
     */
    static class XTree
        extends SuffixTree
    {
        public XTree( int minMatchLen )
        {
            mMinMatchLen = minMatchLen;
        }

        public void decorate()
        {
            mMap = new HashMap<SuffixTreeNode,Set<CharSequence>>();
            walkToDepth(0,mRoot);
        }

        public Set<CharSequence> match( CharSequence probe )
        {
            Set<CharSequence> result = new HashSet<CharSequence>();

            if ( probe.length() >= mMinMatchLen )
            {
                SuffixTreePosition pos = new SuffixTreePosition(mRoot);
                int strLen = probe.length();
                int jjj = 0;
                for ( int iii = 0; iii < strLen; ++iii )
                {
                    if ( jjj < iii )
                        jjj = iii;

                    while ( jjj < strLen && pos.step(probe.charAt(jjj)) )
                        jjj += 1;

                    if ( jjj - iii >= mMinMatchLen )
                    {
                        Set<CharSequence> strings = mMap.get(pos.getNode());
                        if ( strings != null )
                        {
                            result.addAll(strings);
                        }
                    }

                    pos.toSuffix();
                    pos.skip(probe,iii+1,jjj);
                }
            }

            return result;
        }

        private void walkToDepth( int curDepth, SuffixTreeNode node )
        {
            curDepth += node.getLabelLength();
            if ( curDepth >= mMinMatchLen )
            {
                Set<CharSequence> strings = new HashSet<CharSequence>(DEFAULT_INITIAL_CAPACITY);
                buildStringSet(node,strings);
                buildMap(node,strings);
            }
            else if ( !node.isLeaf() )
            {
                for ( char child : CHILDREN )
                {
                    SuffixTreeNode childNode = node.getChild(child);
                    if ( childNode != null )
                    {
                        walkToDepth(curDepth,childNode);
                    }
                }
            }
        }

        private void buildMap( SuffixTreeNode node, Set<CharSequence> strings )
        {
            mMap.put(node,strings);
            if ( !node.isLeaf() )
            {
                for ( char child : CHILDREN )
                {
                    SuffixTreeNode childNode = node.getChild(child);
                    if ( childNode != null )
                    {
                        buildMap(childNode,strings);
                    }
                }
            }
        }

        private static void buildStringSet( SuffixTreeNode node, Set<CharSequence> strings )
        {
            Iterator<CharSequence> itr = node.getStringIterator();
            while ( itr.hasNext() )
            {
                CharSequence str = itr.next();
                if ( str instanceof CharSubSequence )
                {
                    str = ((CharSubSequence)str).getOriginalSequence();
                }
                strings.add(str);
            }
            if ( !node.isLeaf() )
            {
                for ( char child : CHILDREN )
                {
                    SuffixTreeNode childNode = node.getChild(child);
                    if ( childNode != null )
                    {
                        buildStringSet(childNode,strings);
                    }
                }
            }
        }

        private int mMinMatchLen;
        private Map<SuffixTreeNode,Set<CharSequence>> mMap;
        private static final int DEFAULT_INITIAL_CAPACITY = 4;
    }
}
