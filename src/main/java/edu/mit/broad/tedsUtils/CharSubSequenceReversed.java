/*
 * $Id: CharSubSequenceReversed.java 73459 2008-10-10 16:10:03Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

/**
 * A CharSubSequence that presents its characters in reverse order.
 *
 * @author tsharpe
 * @version $Revision: 48897 $
 */
public class CharSubSequenceReversed
    extends CharSubSequence
{
    public CharSubSequenceReversed( CharSequence seq )
    {
        super(seq,0,seq.length());
    }

    public CharSubSequenceReversed( CharSequence seq, int start, int end )
    {
        super(seq,start,end);
    }

    @Override
    public char charAt( int idx )
    {
        return super.charAt(length()-1-idx);
    }

    @Override
    public CharSubSequenceReversed subSequence( int start, int end )
    {
        int originalEnd = getStart() + length();
        return new CharSubSequenceReversed(getOriginalSequence(),originalEnd-end,originalEnd-start);
    }

    @Override
    public String toString()
    {
        int len = length();
        StringBuilder sb = new StringBuilder(len);

        while ( --len >= 0 )
        {
            sb.append(super.charAt(len));
        }

        return sb.toString();
    }
}
