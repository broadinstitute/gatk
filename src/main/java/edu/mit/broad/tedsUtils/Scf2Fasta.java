/**
 * $Id: Scf2Fasta.java 69401 2008-07-03 22:12:03Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */

package edu.mit.broad.tedsUtils;

import java.io.File;
import java.io.IOException;

import edu.mit.broad.tedsUtils.FileTypeFilter.Processor;

/**
 * Reads SCF's and dumps their base calls.
 *
 * @author tsharpe
 * @version $Revision$
 */
public class Scf2Fasta
    implements Processor
{
    public static void main( String[] args )
    {
        new SCFFilter(new Scf2Fasta()).process(args);
    }

    public void process( File scfFile )
    {
        try
        {
            Reading reading = Reading.readSCF(scfFile);
            String seq = reading.getSequence();
            if ( seq == null )
            {
                System.err.println("No sequence: " + scfFile);
            }
            else
            {
                System.out.println(">" + scfFile.getName());
                int len = seq.length();
                int start = 0;
                int end = LINE_LEN;
                while ( end <= len )
                {
                    System.out.println(seq.substring(start,end));
                    start = end;
                    end += LINE_LEN;
                }
                if ( start < len )
                {
                    System.out.println(seq.substring(start,len));
                }
            }
        }
        catch ( IOException e )
        {
            System.err.println("Can't read file: " + scfFile + " -- " + e.getMessage());
        }
    }

    private static final int LINE_LEN = SysProp.iVal("LineLen",60);
}
