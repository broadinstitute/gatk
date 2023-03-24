/*
 * $Id: FastaExtractor.java 69702 2008-07-09 22:07:26Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package edu.mit.broad.tedsUtils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.StringTokenizer;
import java.util.Map.Entry;

/**
 * Make an extract of a FASTA file.
 *
 * @author tsharpe
 * @version $Revision: 69702 $
 */
public class FastaExtractor
{
    public static void main( String[] args )
        throws IOException
    {
        FastaExtractor fe = new FastaExtractor();
        BufferedOutputStream bos = new BufferedOutputStream(System.out);
        BufferedReader regionsRdr = new BufferedReader(new InputStreamReader(System.in));
        String line;
        while ( (line = regionsRdr.readLine()) != null )
        {
            StringTokenizer tok = new StringTokenizer(line," \t");
            String chrName = tok.nextToken().intern();
            int begin = Integer.parseInt(tok.nextToken()) - 1;
            int end = Integer.parseInt(tok.nextToken());
            String comment;
            if ( tok.hasMoreTokens() )
            {
                comment = tok.nextToken();
            }
            else
            {
                comment = chrName + ":" + (begin+1) + "-" + end;
            }
            fe.extract( bos, comment, chrName, begin, end );
        }

        bos.close();
    }

    public FastaExtractor()
        throws IOException
    {
        this(PATHS_FILE);
    }

    public FastaExtractor( String pathsFile )
        throws IOException
    {
        this(pathsFile,DEFAULT_MAX_FILES_OPEN);
    }

    public FastaExtractor( String pathsFile, int maxFilesOpen )
        throws IOException
    {
        BufferedReader pathsRdr = new BufferedReader(new FileReader(pathsFile));
        HashMap<String,FastaFileExtractor> map = new HashMap<String,FastaFileExtractor>();
        String line;
        while ( (line = pathsRdr.readLine()) != null )
        {
            StringTokenizer tok = new StringTokenizer(line);
            String chrName = tok.nextToken().intern();
            String path = tok.nextToken();
            map.put(chrName,new FastaFileExtractor(path));
        }
        pathsRdr.close();
        mFileExtractorMap = map;
        mCacheMap = new CacheMap(maxFilesOpen);
    }

    public void extract( BufferedOutputStream bos, String comment, String chrName, int begin, int end )
        throws IOException
    {
        FastaFileExtractor ffe = getFFE(chrName);
        if ( ffe == null )
        {
            System.err.println("Can't find FASTA file for " + chrName);
        }
        else
        {
            if ( comment != null )
            {
                ffe.comment(bos,comment);
            }
            ffe.extract(bos,begin,end);
        }
    }

    public String getSequence( String chrName, int begin, int end )
        throws IOException
    {
        FastaFileExtractor ffe = getFFE(chrName);
        if ( ffe == null )
        {
            throw new IOException("Can't find FASTA file for " + chrName);
        }

        return ffe.getSequence(begin,end);
    }

    private FastaFileExtractor getFFE( String chrName )
    {
        FastaFileExtractor result = mFileExtractorMap.get(chrName);
        if ( result != null )
        {
            // manage the number of open files
            mCacheMap.put(chrName,result);
        }
        return result;
    }

    private HashMap<String,FastaFileExtractor> mFileExtractorMap;
    private CacheMap mCacheMap;

    private static final String PATHS_FILE = SysProp.sVal("PathsFile","paths");
    private static final int DEFAULT_MAX_FILES_OPEN = SysProp.iVal("MaxFilesOpen",3);

    public static class FastaFileExtractor
    {
        public FastaFileExtractor( String path )
        {
            mPath = new File(path);
        }

        public String getSequence( int begin, int end )
            throws IOException
        {
            initFile();
            int pos1 = mInitialOffset + begin / mLineLength + begin;
            int pos2 = mInitialOffset + end / mLineLength + end;
            int len = pos2 - pos1;
            int totalLineLen = mLineLength + 1;
            int lfModulus = mLineLength - (pos1 - mInitialOffset) % (totalLineLen);

            byte[] buf = new byte[len];
            mFile.seek(pos1);
            mFile.readFully(buf);

            StringBuilder sb = new StringBuilder();
            for ( int idx = 0; idx < len; ++idx )
            {
                if ( idx % totalLineLen != lfModulus )
                {
                    sb.append((char)buf[idx]);
                }
            }

            return sb.toString();
        }

        public int getMaxEnd()
        {
            return (int)(mLineLength * (mPath.length() - mInitialOffset))/(mLineLength + 1);
        }

        public void comment( BufferedOutputStream bos, String comment )
            throws IOException
        {
            bos.write('>');
            bos.write(comment.getBytes());
            bos.write('\n');
        }

        public void extract( BufferedOutputStream bos, int begin, int end )
            throws IOException
        {
            initFile();
            int coreBegin = 0;
            int coreEnd = end - begin;
            if ( FLANKING_SEQUENCE > 0 )
            {
                coreBegin = Math.min(begin,FLANKING_SEQUENCE);
                begin -= coreBegin;
                coreEnd = end - begin;
                end = Math.min(end+FLANKING_SEQUENCE,getMaxEnd());
            }
            mFile.seek(mInitialOffset + begin / mLineLength + begin );

            byte[] buf = new byte[8192];
            int col = 0;
            int len = end - begin;
            int cnt = 0;
            while ( len > 0 )
            {
                int bufLen = mFile.read(buf);
                if ( bufLen == -1 )
                {
                    throw new IOException("Premature EOF");
                }

                for ( int idx = 0; idx < bufLen; ++idx )
                {
                    int bbb;
                    if ( (bbb = buf[idx]) != '\n' )
                    {
                        if ( cnt < coreBegin || cnt >= coreEnd )
                        {
                            bbb = Character.toLowerCase(bbb);
                        }
                        else if ( UPCASE_EXTRACT )
                        {
                            bbb = Character.toUpperCase(bbb);
                        }
                        cnt += 1;
                        bos.write(bbb);
                        if ( ++col % OUTPUT_LINE_LEN == 0 )
                        {
                            bos.write('\n');
                        }
                        if ( --len == 0 )
                        {
                            break;
                        }
                    }
                }
            }

            if ( col % OUTPUT_LINE_LEN != 0 )
            {
                bos.write('\n');
            }
        }

        public void getSkinny()
        {
            if ( mFile != null )
            {
                try
                {
                    mFile.close();
                }
                catch ( IOException ioe )
                {
                    System.err.println("Unable to close file: " + ioe.getMessage());
                }
                mFile = null;
            }
        }

        private void initFile()
            throws IOException
        {
            if ( mLineLength == 0 )
            {
                init();
            }

            if ( mFile == null )
            {
                mFile = new RandomAccessFile(mPath,"r");
            }
        }

        private void init()
            throws IOException
        {
            BufferedInputStream bis = new BufferedInputStream(new FileInputStream(mPath));
            int bbb;
            while ( (bbb = bis.read()) == '>' )
            {
                ++mInitialOffset;

                while ( (bbb = bis.read()) != -1 )
                {
                    ++mInitialOffset;
                    if ( bbb == '\n' )
                        break;
                }
            }
            while ( bbb != -1 )
            {
                if ( bbb == '\n' )
                    break;
                ++mLineLength;
                bbb = bis.read();
            }
        }

        private File mPath;
        private RandomAccessFile mFile;
        private int mLineLength;
        private int mInitialOffset;

        private static final int FLANKING_SEQUENCE = SysProp.iVal("FlankingSequence",0);
        private static final int OUTPUT_LINE_LEN = SysProp.iVal("LineLength",60);
        private static final boolean UPCASE_EXTRACT = SysProp.bVal("Upcase",false);
    }

    static class CacheMap
        extends LinkedHashMap<String,FastaFileExtractor>
    {
        public CacheMap( int nEntries )
        {
            super(nEntries+1);
            mNEntries = nEntries;
        }

        @Override
        protected boolean removeEldestEntry( Entry<String,FastaFileExtractor> entry )
        {
            boolean result = false;
            if ( size() > mNEntries )
            {
                entry.getValue().getSkinny();
                result = true;
            }
            return result;
        }

        private int mNEntries;
        private static final long serialVersionUID = 630916669721404447L;
    }
}
