package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.*;
import org.apache.hadoop.io.NullWritable;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.seqdoop.hadoop_bam.BAMRecordWriter;
import org.seqdoop.hadoop_bam.SAMRecordWritable;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;

@CommandLineProgramProperties(summary="Gather clustered split reads using spark",
        oneLineSummary="Gather clustered split reads using spark",
        programGroup = SparkProgramGroup.class)
public class GatherSplitReadsSpark extends GATKSparkTool
{
    private static final long serialVersionUID = 1L;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    private static final class EventEvidence implements Comparable<EventEvidence>
    {
        public EventEvidence( final int someLocus )
        {
            locus = someLocus;
            counter = 0;
        }

        public EventEvidence( final int someLocus, final long someCounter )
        {
            locus = someLocus;
            counter = someCounter;
        }

        public int getLocus() { return locus; }

        public int compareTo( final EventEvidence ee )
        {
            int result = Integer.compare(locus,ee.locus);
            if ( result == 0 ) result = Long.compare(counter,ee.counter);
            return result;
        }

        private final int locus;
        private final long counter; // used as a uniquifier
    }

    private static final class SplitReadClusterer implements Iterator<GATKRead>
    {
        public SplitReadClusterer( final Iterator<GATKRead> someInputIterator )
        {
            inputIterator = someInputIterator;
        }

        @Override
        public boolean hasNext()
        {
            while ( !outputIterator.hasNext() )
            {
                if ( !inputIterator.hasNext() )
                    return false;
                outputIterator = processRead(inputIterator.next());
            }
            return true;
        }

        @Override
        public GATKRead next()
        {
            if ( !hasNext() )
                throw new NoSuchElementException("Iterator<GATKRead> is exhausted.");
            return outputIterator.next();
        }

        private Iterator<GATKRead> processRead( final GATKRead read )
        {
            if ( !read.failsVendorQualityCheck() && !read.isUnmapped() &&
                    !read.isDuplicate() && read.getMappingQuality() > 0 &&
                    read.getStart() != ReadConstants.UNSET_POSITION )
            {
                if ( partitionContig == null )
                    partitionContig = read.getContig();
                else if ( partitionContig != read.getContig() )
                    throw new GATKException("This tool cannot handle partitions that cross contig boundaries.  "+
                            "Please reindex your input to provide single-contig partitions.");

                final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
                int matchLen = 0;
                for ( CigarElement ele : cigarElements )
                    if (ele.getOperator() == CigarOperator.MATCH_OR_MISMATCH)
                        matchLen += ele.getLength();
                if ( matchLen >= MIN_MATCH_LEN )
                {
                    final ListIterator<CigarElement> itr0 = cigarElements.listIterator();
                    if ( itr0.hasNext() )
                    {
                        CigarElement firstEle = itr0.next();
                        if ( firstEle.getOperator() == CigarOperator.HARD_CLIP && itr0.hasNext() )
                            firstEle = itr0.next();
                        if ( firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                                firstEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                                highQualityRegion(read.getBaseQualities(), 0) )
                            return cluster(read, read.getStart()); // NON-STRUCTURED return
                    }
                    final ListIterator<CigarElement> itrN = cigarElements.listIterator(cigarElements.size());
                    if ( itrN.hasPrevious() )
                    {
                        CigarElement lastEle = itrN.previous();
                        if ( lastEle.getOperator() == CigarOperator.HARD_CLIP && itrN.hasPrevious() )
                            lastEle = itrN.previous();
                        if ( lastEle.getOperator() == CigarOperator.SOFT_CLIP &&
                                lastEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                                highQualityRegion(read.getBaseQualities(), read.getLength()-lastEle.getLength()) )
                            return cluster(read, read.getEnd()); // NON-STRUCTURED return
                    }
                }
            }
            return Collections.emptyIterator();
        }

        private Iterator<GATKRead> cluster( final GATKRead read, final int locus )
        {
            final EventEvidence newEvidence = new EventEvidence(locus,++counter);
            locMap.put(newEvidence,read);

            // clean out old stuff that can't possibly be interesting anymore
            final Iterator<Map.Entry<EventEvidence,GATKRead>> itr = locMap.entrySet().iterator();
            while ( itr.hasNext() )
            {
                final EventEvidence oldEvidence = itr.next().getKey();
                if ( oldEvidence.getLocus() + MAX_LOCUS_DIST > newEvidence.getLocus() )
                    break;
                itr.remove();
            }

            // find all the evidence in a window surrounding the locus of the new evidence
            final SortedMap<EventEvidence,GATKRead> windowMap = locMap
                    .tailMap(new EventEvidence(newEvidence.getLocus()-CLUSTER_WINDOW))
                    .headMap(new EventEvidence(newEvidence.getLocus()+CLUSTER_WINDOW + 1));
            if ( windowMap.size() >= MIN_EVIDENCE )
            {
                final List<GATKRead> list = new LinkedList<>();
                for ( Map.Entry<EventEvidence,GATKRead> entry : windowMap.entrySet() )
                {
                    if (entry.getValue() != null)
                    {
                        list.add(entry.getValue());
                        entry.setValue(null);
                    }
                }
                return list.iterator(); // NON-STRUCTURED return
            }
            return Collections.emptyIterator();
        }

        private static boolean highQualityRegion( final byte[] quals, int idx )
        {
            for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx )
                if ( quals[idx] < MIN_QUALITY )
                    return false; // NON-STRUCTURED return
            return true;
        }

        private String partitionContig = null;
        private long counter = 0;
        private final Iterator<GATKRead> inputIterator;
        private final SortedMap<EventEvidence,GATKRead> locMap = new TreeMap<>();
        private Iterator<GATKRead> outputIterator = Collections.emptyIterator();

        private static final int MIN_MATCH_LEN = 45; // minimum length of matched portion of an interesting alignment
        private static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
        private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window
        private static final int MAX_LOCUS_DIST = 500; // stale evidence distance, should be somewhat longer than a read
        private static final int CLUSTER_WINDOW = 2; // size of locus window in which we cluster event evidence
        private static final int MIN_EVIDENCE = 3; // minimum evidence count in a cluster
    }

    private static final class SAMWindowSorter implements Iterator<GATKRead>, Iterable<GATKRead>
    {
        public SAMWindowSorter( final Iterator<GATKRead> someInputIterator, final int someWindowSize )
        {
            inputIterator = someInputIterator;
            windowSize = someWindowSize;
            fillBuffer();
        }

        @Override
        public boolean hasNext()
        {
            return !recordSet.isEmpty();
        }

        @Override
        public GATKRead next()
        {
            final Iterator<GATKRead> itr = recordSet.iterator();
            final GATKRead result = itr.next();
            itr.remove();
            fillBuffer();
            return result;
        }

        @Override
        public Iterator<GATKRead> iterator() { return this; }

        private void fillBuffer()
        {
            while ( getSetSpan() <= windowSize && inputIterator.hasNext() )
                recordSet.add(inputIterator.next());
        }

        private int getSetSpan()
        {
            if ( recordSet.isEmpty() )
                return 0;
            final GATKRead first = recordSet.first();
            final GATKRead last = recordSet.last();
            return last.getStart() - first.getStart();
        }

        private static final class ReadComparator implements Comparator<GATKRead>
        {
            @Override
            public int compare( final GATKRead rec1, final GATKRead rec2 )
            {
                int result = Integer.compare(rec1.getStart(), rec2.getStart());
                if ( result == 0 ) result = rec1.getName().compareTo(rec2.getName());
                if ( result == 0 ) result = -Boolean.compare(rec1.isFirstOfPair(),rec2.isFirstOfPair());
                return result;
            }
        }

        private final SortedSet<GATKRead> recordSet = new TreeSet<>(new ReadComparator());
        private final Iterator<GATKRead> inputIterator;
        private final int windowSize;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx )
    {
        final Broadcast<SAMFileHeader> header = ctx.broadcast(getHeaderForReads());
        if ( header.value().getSortOrder() != SAMFileHeader.SortOrder.coordinate )
            throw new IllegalStateException("The BAM must be coordinate sorted.");

        final JavaRDD<GATKRead> clusteredReads =
            getReads().mapPartitions(readItr ->
                { return new SAMWindowSorter(new SplitReadClusterer(readItr),SAM_WINDOW_SIZE); }, true);

        try
        { ReadsSparkSink.writeReads(ctx, output, clusteredReads, header.value(), ReadsWriteFormat.SINGLE); }
        catch ( IOException e )
        { throw new GATKException("Unable to write BAM" + output, e); }
    }

    private static final int SAM_WINDOW_SIZE = 500;
}
