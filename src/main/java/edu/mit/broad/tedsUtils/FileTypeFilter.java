/*
 * $Id: FileTypeFilter.java 70621 2008-07-29 18:51:48Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

import java.io.File;
import java.io.FileFilter;
import java.util.regex.Pattern;

/**
 * Processes a list of file and directory names, looking for files whose names match a pattern.
 * If you specify -DRecursiveDirProcessing=true, directories are scanned recursively, otherwise,
 * and by default, only the named directory itself is scanned for files matching the pattern.
 * @author tsharpe
 */
public class FileTypeFilter
    implements FileFilter
{
    public interface Processor
    {
        void process( File file );
    }

    public FileTypeFilter( Pattern pattern, Processor processor )
    {
        this(pattern,processor,RECURSIVE_DIR_PROCESSING);
    }

    public FileTypeFilter( Pattern pattern, Processor processor, boolean recursiveDirProcessing )
    {
        mPattern = pattern;
        mProcessor = processor;
        mRecursiveDirProcessing = recursiveDirProcessing;
    }

    public void process( String[] args )
    {
        for ( String arg : args )
        {
            process(arg);
        }
    }

    public void process( String fileName )
    {
        File file = new File(fileName);
        if ( file.isDirectory() )
        {
            file.listFiles(this);
        }
        else if ( !file.exists() )
        {
            System.err.println("No such file or directory: " + file.getPath());
        }
        else
        {
            accept(file);
        }
    }

    public boolean accept( File file )
    {
        if ( mPattern.matcher(file.getName()).matches())
        {
            mProcessor.process(file);
        }
        else if ( mRecursiveDirProcessing && file.isDirectory() )
        {
            file.listFiles(this);
        }
        return false;
    }

    private Pattern mPattern;
    private Processor mProcessor;
    private boolean mRecursiveDirProcessing;

    private static final boolean RECURSIVE_DIR_PROCESSING = SysProp.bVal("RecursiveDirProcessing",false);
}
