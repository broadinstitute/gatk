package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
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

import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(summary="Gather clustered split reads using spark",
        oneLineSummary="Gather clustered split reads using spark",
        programGroup = SparkProgramGroup.class)
public class GatherSplitReadsSpark extends GATKSparkTool
{
    private static final long serialVersionUID = 1L;
    private static final int OUTPUT_SORT_WINDOW_SIZE = 500;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Override
    public boolean requiresReads()
    {
        return true;
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
                new WindowSorter(new BreakpointClusterer(new SplitReadDetector(),
                                                         new DiscordantPairDetector(),
                                                         readItr),OUTPUT_SORT_WINDOW_SIZE), true)
            .coalesce(1)
            .mapPartitions(readItr ->
                new WindowSorter(readItr,OUTPUT_SORT_WINDOW_SIZE), true);

        try
        { ReadsSparkSink.writeReads(ctx,output,clusteredReads,header,ReadsWriteFormat.SINGLE); }
        catch ( IOException e )
        { throw new GATKException("Unable to write BAM" + output, e); }
    }

    private static final class EventLocus implements Comparable<EventLocus>
    {
        public EventLocus( final int someLocusStart, final int intervalWidth )
        {
            locusStart = someLocusStart;
            locusEnd = someLocusStart + intervalWidth;
        }

        public int getLocusStart() { return locusStart; }
        public int getLocusEnd() { return locusEnd; }

        public EventLocus makeUnique( long someId )
        {
            if ( id == 0 )
                id = someId;
            return this;
        }

        public int compareTo( final EventLocus ee )
        {
            int result = Integer.compare(locusStart,ee.locusStart);
            if ( result == 0 ) result = Integer.compare(locusEnd,ee.locusEnd);
            if ( result == 0 ) result = Long.compare(id,ee.id);
            return result;
        }

        private final int locusStart;
        private final int locusEnd;
        private long id = 0; // used as a uniquifier
    }

    private static final class SplitReadDetector
    {
        public EventLocus isSplit( GATKRead read )
        {
            final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
            final ListIterator<CigarElement> itr0 = cigarElements.listIterator();
            if ( itr0.hasNext() )
            {
                CigarElement firstEle = itr0.next();
                if ( firstEle.getOperator() == CigarOperator.HARD_CLIP && itr0.hasNext() )
                    firstEle = itr0.next();
                if ( firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                        firstEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                        highQualityRegion(read.getBaseQualities(), 0) )
                    return new EventLocus(read.getStart()-SOFT_CLIP_LOCUS_WIDTH/2, SOFT_CLIP_LOCUS_WIDTH);
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
                    return new EventLocus(read.getEnd()-SOFT_CLIP_LOCUS_WIDTH/2, SOFT_CLIP_LOCUS_WIDTH);
            }

            return null;
        }

        private static boolean highQualityRegion( final byte[] quals, int idx )
        {
            for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx )
                if ( quals[idx] < MIN_QUALITY )
                    return false;

            return true;
        }

        private static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
        private static final int SOFT_CLIP_LOCUS_WIDTH = 4; // uncertainty in event locus for soft clip
        private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window
    }

    private static final class DiscordantPairDetector
    {
        public EventLocus isDiscordant( GATKRead read )
        {
            if ( !read.mateIsUnmapped() && read.getContig() == read.getMateContig() )
            {
                if ( read.isReverseStrand() == read.mateIsReverseStrand() || !read.isProperlyPaired() )
                {
                    int locus = read.getFragmentLength() < 0 ?
                            read.getStart()-FUNKY_PAIR_LOCUS_WIDTH :
                            read.getEnd();
                    return new EventLocus(locus, FUNKY_PAIR_LOCUS_WIDTH);
                }
            }

            return null;
        }

        private static final int FUNKY_PAIR_LOCUS_WIDTH = 100; // uncertainty in event locus for discordant pair
    }

    private static final class BreakpointClusterer implements Iterator<GATKRead>
    {
        public BreakpointClusterer( final SplitReadDetector someSplitReadDetector,
                                    final DiscordantPairDetector someDiscordantPairDetector,
                                    final Iterator<GATKRead> someInputIterator )
        {
            splitReadDetector = someSplitReadDetector;
            discordantPairDetector = someDiscordantPairDetector;
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
                    EventLocus eventLocus = splitReadDetector.isSplit(read);
                    if ( eventLocus != null )
                        return cluster(eventLocus,read);
                    eventLocus = discordantPairDetector.isDiscordant(read);
                    if ( eventLocus != null )
                        return cluster(eventLocus,read);
                }
            }
            return Collections.emptyIterator();
        }

        private Iterator<GATKRead> cluster( final EventLocus newEventLocus, final GATKRead read )
        {
            final int locusStart = newEventLocus.getLocusStart();
            final int locusEnd = newEventLocus.getLocusEnd();

            // clean out old stuff that can't possibly be interesting anymore
            final Iterator<Map.Entry<EventLocus,GATKRead>> itr = locMap.entrySet().iterator();
            final int staleEnd = locusStart - MAX_LOCUS_DIST;
            while ( itr.hasNext() )
            {
                final EventLocus eventLocus = itr.next().getKey();
                if ( eventLocus.getLocusStart() >= locusStart )
                    break;
                if ( eventLocus.getLocusEnd() <= staleEnd )
                    itr.remove();
            }

            locMap.put(newEventLocus.makeUnique(++eventIdCounter),read);

            // see if there are a sufficient number of overlapping putative events to call a cluster
            final List<Map.Entry<EventLocus,GATKRead>> overlappers = new LinkedList<>();
            for ( Map.Entry<EventLocus,GATKRead> entry : locMap.entrySet() )
            {
                final EventLocus eventLocus = entry.getKey();
                if ( eventLocus.getLocusStart() >= locusEnd )
                    break;
                if ( eventLocus.getLocusEnd() > locusStart )
                    overlappers.add(entry);
            }
            if ( overlappers.size() >= MIN_EVIDENCE )
            {
                final List<GATKRead> list = new LinkedList<>();
                for ( Map.Entry<EventLocus,GATKRead> entry : overlappers )
                {
                    if (entry.getValue() != null)
                    {
                        list.add(entry.getValue());
                        entry.setValue(null);
                    }
                }
                return list.iterator();
            }
            return Collections.emptyIterator();
        }

        private SplitReadDetector splitReadDetector;
        private DiscordantPairDetector discordantPairDetector;
        private final Iterator<GATKRead> inputIterator;
        private String currentContig = null;
        private long eventIdCounter = 0;
        private final SortedMap<EventLocus,GATKRead> locMap = new TreeMap<>();
        private Iterator<GATKRead> outputIterator = Collections.emptyIterator();

        private static final int MIN_MATCH_LEN = 45; // minimum length of matched portion of an interesting alignment
        private static final int MAX_LOCUS_DIST = 300; // stale evidence distance, should be about the fragment length
        private static final int MIN_EVIDENCE = 5; // minimum evidence count in a cluster
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
                return 47*read.getName().hashCode() ^ Boolean.hashCode(isFirst);
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
}
