package htsjdk.samtools;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by tsharpe on 11/25/15.
 */
public class DumpBAMFileIndex {
    public static void main( String[] args )
    {
        try
        {
            BAMFileReader reader = new BAMFileReader(new File(args[0]), new File(args[1]), false, ValidationStringency.SILENT, new DefaultSAMRecordFactory());
            reader.enableIndexCaching(true);
            CachingBAMFileIndex bamIndex = (CachingBAMFileIndex)reader.getIndex();
            int nRefs = reader.getFileHeader().getSequenceDictionary().size();
            for ( int refIdx = 0; refIdx != nRefs; ++refIdx )
            {
                SortedSet<Long> locs = new TreeSet<Long>();
                BAMIndexContent indexContent = bamIndex.getQueryResults(refIdx);
                for ( Chunk chunk : indexContent.getAllChunks() )
                {
                    locs.add(chunk.getChunkStart());
                    locs.add(chunk.getChunkEnd());
                }
                LinearIndex linearIndex = indexContent.getLinearIndex();
                for ( Long loc : linearIndex.getIndexEntries() )
                    locs.add(loc);
                System.out.print("Reference "+refIdx+":\n");
                long lastOff = 0;
                for ( Long loc : locs )
                {
                    long off = loc >> 16;
                    System.out.print(off+"\t"+(loc&0xffff)+"\t"+(off-lastOff)+"\n");
                    lastOff = off;
                }
            }
        }
        catch ( IOException e )
        {
            throw new IllegalStateException("Failed.",e);
        }
    }
}
