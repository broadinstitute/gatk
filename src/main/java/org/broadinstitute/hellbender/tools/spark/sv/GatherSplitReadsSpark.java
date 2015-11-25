package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
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
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;
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
                String readContig = read.getContig();
                if ( currentContig != readContig )
                {
                    currentContig = readContig;
                    locMap.clear();
                }

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

        private String currentContig = null;
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

    private static final class WindowSorter implements Iterator<GATKRead>, Iterable<GATKRead>
    {
        private static final class GATKKey implements Comparable<GATKKey>
        {
            GATKKey( int someRelativeContigNo, GATKRead someRead )
            {
                relativeContigNo = someRelativeContigNo;
                locus = someRead.getStart();
                isFirst = someRead.isFirstOfPair();
                read = someRead;
            }

            @Override
            public int compareTo( GATKKey that )
            {
                int result = Integer.compare(this.relativeContigNo,that.relativeContigNo);
                if ( result == 0 ) result = Integer.compare(this.locus,that.locus);
                if ( result == 0 ) result = this.read.getName().compareTo(that.read.getName());
                if ( result == 0 ) result = -Boolean.compare(this.isFirst,that.isFirst);
                return result;
            }

            @Override
            public boolean equals( Object obj )
            {
                if ( !(obj instanceof GATKKey) ) return false;
                GATKKey that = (GATKKey)obj;
                return this.read.getName() == that.read.getName() && this.isFirst == that.isFirst;
            }

            @Override
            public int hashCode()
            {
                return 47*read.hashCode() ^ Boolean.hashCode(isFirst);
            }

            public int distanceTo( GATKKey that )
            {
                if ( this.relativeContigNo != that.relativeContigNo )
                    return Integer.MAX_VALUE;
                return that.locus - this.locus;
            }

            public GATKRead getRead() { return read; }

            private final int relativeContigNo;
            private final int locus;
            private final boolean isFirst;
            private final GATKRead read;
        }

        public WindowSorter( final Iterator<GATKRead> someInputIterator,
                             final int someWindowSize )
        {
            currentContig = null;
            contigNo = 0;
            recordSet = new TreeSet<>();
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
            final Iterator<GATKKey> itr = recordSet.iterator();
            final GATKRead result = itr.next().getRead();
            itr.remove();
            fillBuffer();
            return result;
        }

        @Override
        public Iterator<GATKRead> iterator() { return this; }

        private void fillBuffer()
        {
            while ( getSetSpan() <= windowSize && inputIterator.hasNext() )
            {
                GATKRead read = inputIterator.next();
                String readContig = read.getContig();
                if ( currentContig == null || currentContig != readContig )
                {
                    currentContig = readContig;
                    contigNo += 1;
                }
                recordSet.add(new GATKKey(contigNo,read));
            }
        }

        private int getSetSpan()
        {
            if ( recordSet.isEmpty() ) return 0;
            final GATKKey first = recordSet.first();
            final GATKKey last = recordSet.last();
            if ( first == last ) return 0;
            return first.distanceTo(last);
        }

        private String currentContig;
        private int contigNo;
        private final SortedSet<GATKKey> recordSet;
        private final Iterator<GATKRead> inputIterator;
        private final int windowSize;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx )
    {
        final SAMFileHeader header = getHeaderForReads();
        if ( header.getSortOrder() != SAMFileHeader.SortOrder.coordinate )
            throw new GATKException("The BAM must be coordinate sorted.");

        final JavaRDD<GATKRead> clusteredReads =
            getReads()
                .mapPartitions(readItr ->
                { return new WindowSorter(new SplitReadClusterer(readItr),SAM_WINDOW_SIZE); }, true)
                .coalesce(1)
                .mapPartitions(readItr ->
                { return new WindowSorter(readItr,SAM_WINDOW_SIZE); }, true);

        try
        { ReadsSparkSink.writeReads(ctx,output,clusteredReads,header,ReadsWriteFormat.SINGLE); }
        catch ( IOException e )
        { throw new GATKException("Unable to write BAM" + output, e); }
    }

    private static final int SAM_WINDOW_SIZE = 500;
}
