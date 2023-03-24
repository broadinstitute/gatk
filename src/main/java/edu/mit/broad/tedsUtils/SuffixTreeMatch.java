/*
 * $Id: SuffixTreeMatch.java 73459 2008-10-10 16:10:03Z tsharpe $
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

/**
 * A part of a probe and a part of a target sequence that are identical.
 *
 * @author tsharpe
 * @version $Revision: 13430 $
 */
public class SuffixTreeMatch
{
    public SuffixTreeMatch( CharSubSequence probe, CharSubSequence target )
    {
        mProbe = probe;
        mTarget = target;
    }

    public CharSubSequence getProbe()
    {
        return mProbe;
    }
    
    public CharSubSequence getTarget()
    {
        return mTarget;
    }

    private CharSubSequence mProbe;
    private CharSubSequence mTarget;
}
