/*
 * $Id: RefReader.java 40304 2007-06-07 15:34:19Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Reads fasta into a String.
 * Ignores comments lines.
 * Very forgiving:  just ignores characters that can't be interpreted as nucleotides.
 * @author tsharpe
 * @version $Revision: 40304 $
 */
public class RefReader
{
    /**
     * Reads from stdin.
     */
    public RefReader()
        throws IOException
    {
        this(null);
    }

    /**
     * Reads from file.
     * @param fileName File name.
     * @throws IOException
     */
    public RefReader( String fileName )
        throws IOException
    {
        BufferedReader rdr = null;
        try
        {
            rdr = new BufferedReader(fileName==null ? new InputStreamReader(System.in) : new FileReader(fileName));
            init(rdr);
        }
        finally
        {
            if ( rdr != null )
            {
                rdr.close();
            }
        }
    }

    /**
     * Returns the reference sequence.
     */
    public String getRef()
    {
        return mRef;
    }

    /**
     * Returns the reverse-complemented reference sequence.
     */
    public String getRevRef()
    {
        if ( mRevRef == null && mRef != null )
        {
            mRevRef = FASTAUtils.reverseComplement(mRef);
        }
        return mRevRef;
    }

    /**
     * Returns the last comment line with the '>' trimmed off.
     */
    public String getComment()
    {
        return mComment;
    }

    private void init( BufferedReader rdr )
        throws IOException
    {
        StringBuilder sb = new StringBuilder();
        String line;
        while ( (line = rdr.readLine()) != null )
        {
            int len = line.length();
            if ( len > 0 )
            {
                if ( line.charAt(0) == '>' )
                {
                    mComment = line.substring(1).trim();
                }
                else
                {
                    for ( int idx = 0; idx < len; ++idx )
                    {
                        char chr = line.charAt(idx);
                        if ( "ACGTNactgn".indexOf(chr) != -1 )
                        {
                            sb.append(chr);
                        }
                    }
                }
            }
        }
        mRef = sb.toString();
    }

    private String mRef;
    private String mRevRef;
    private String mComment;
}
