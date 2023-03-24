/*
 * $Id: SCFFilter.java 43469 2007-08-23 21:57:42Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

import java.util.regex.Pattern;

/**
 * Handles processing of a list of SCF files.
 * Also processes directories of SCFs, and can do so recursively.
 * On discovering an SCF, it hands off the actual processing to an
 * implementation of FileTypeFilter.Processor.
 *
 * @author tsharpe
 * @version $Revision: 43469 $
 */
public class SCFFilter
    extends FileTypeFilter
{
    public SCFFilter( Processor processor )
    {
        super(gPattern,processor);
    }

    public SCFFilter( Processor processor, boolean recursiveDirProcessing )
    {
        super(gPattern,processor,recursiveDirProcessing);
    }

    private static final Pattern gPattern = Pattern.compile(".*[.]scf(?:[.]gz)?+");
}
