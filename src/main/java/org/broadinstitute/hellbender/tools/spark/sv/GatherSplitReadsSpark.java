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
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import scala.Tuple2;

import java.io.*;
import java.util.*;

@CommandLineProgramProperties(summary="Produce FASTQs of breakpoint candidates",
        oneLineSummary="Produce mini-FASTQs of reads sharing kmers with putative breakpoints for local assembly",
        programGroup = SparkProgramGroup.class)
public class GatherSplitReadsSpark extends GATKSparkTool
{
    private static final long serialVersionUID = 1L;
    private static final int KLEN = 63;

    @Argument(doc = "dir for output FASTQs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputDir;

    // This is a path to a file of kmers that appear too frequently in the reference to be usable as probes to localize
    // reads.  We don't calculate it here, because it depends only on the reference.
    // The program FindBadGenomicKmers can produce such a list for you.
    @Argument(doc = "ubiquitous kmer list", shortName = "KS", fullName = "kmerKillSet", optional = false)
    private String kmerKillSet;

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

        // Process the input data to produce a map of Kmer -> List of breakpoint IDs.
        // (A breakpoint ID is just an arbitrary identifier for the putative breakpoint.)
        // See the doc for findClusters to see a description of how this big Spark pipeline works.
        final Broadcast<Map<Kmer,List<Long>>> broadcastMap = ctx.broadcast(findClusters(ctx));

        // Process the input again, this time pulling all reads that contain kmers associated with a given breakpoint,
        // and writing those reads into a separate FASTQ for each breakpoint.
        getGoodReads()
            .mapPartitionsToPair(readItr -> new ReadPuller(readItr,broadcastMap.value()))
            .groupByKey()
            .foreach(tupl -> writeFastq(tupl._1,tupl._2.iterator(),outputDir));
    }

    public JavaRDD<GATKRead> getGoodReads() {
        return getUnfilteredReads().filter(read ->
            !read.isSecondaryAlignment() && !read.isSupplementaryAlignment() &&
                    !read.isDuplicate() && !read.failsVendorQualityCheck());
    }

    // Read a file of kmers to ignore.
    // Each line must be exactly KLEN characters long, and must match [ACGT]*.
    private Set<Kmer> readKillSet()
    {
        final Set<Kmer> killSet = new HashSet<>();

        try
        {
            final BufferedReader rdr = new BufferedReader(new FileReader(kmerKillSet));
            String line;
            while ( (line = rdr.readLine()) != null )
            {
                if ( line.length() != KLEN )
                    throw new GATKException("Kmer kill set contains a line of length " + line.length() +
                            " but we were expecting K=" + KLEN);

                final StringKmerizer kmerizer = new StringKmerizer(line, KLEN);
                if ( !kmerizer.hasNext() )
                    throw new GATKException("Unable to kmerize the kmer kill set string '" + line + "'.");

                killSet.add(kmerizer.next());
            }
        }
        catch ( IOException e )
        {
            throw new GATKException("Unable to read kmer kill set.", e);
        }

        return killSet;
    }

    // Produce a map of Kmer -> List of unique breakpoint IDs.
    // This is done by looking for reads that appear to be split, or for reads that are part of discordant pairs.
    // When a sufficient number of these reads are mapped onto reference near one another, those reads are taken
    // together as evidence of a putative breakpoint.
    // Each such breakpoint is numbered, and the reads that supported that breakpoint are kmerized.
    // Kmers that appear on a kill-list are ignored.
    // This produces a map of Breakpoint ID -> List of Kmers, which is inverted to produce the final output.
    private Map<Kmer,List<Long>> findClusters( final JavaSparkContext ctx )
    {
        // read the ignored kmer list
        final Broadcast<Set<Kmer>> killSet = ctx.broadcast(readKillSet());

        // produce a set of Breakpoint IDs for each of an interesting set of Kmers by...
        // 1) identifying putative breakpoints supported by a cluster of funky reads
        // 2) kmerizing all the funky reads in the cluster
        // 3) grouping by Kmer
        final List<Tuple2<Kmer,Set<Long>>> readPullingKmers =
            getGoodReads()
                .mapPartitionsWithIndex((partIdx,readItr) ->
                    new ClusterKmerizer(
                        partIdx,
                        new BreakpointClusterer(
                                new SplitReadDetector(),
                                new DiscordantPairDetector(),
                                readItr),
                        killSet.value()),
                    false)
                .mapToPair(tupl -> tupl)
                .combineByKey(breakId -> { Set<Long> coll = new HashSet<>(4); coll.add(breakId); return coll; },
                        (coll,breakId) -> { coll.add(breakId); return coll; },
                        (coll1,coll2) -> { coll1.addAll(coll2); return coll1; })
                .collect();

        // turn the list of interesting Kmers and the set of putative breakpoints to which they belong into a map
        int nKmers = readPullingKmers.size();
        final Map<Kmer,List<Long>> readPullingMap = new HashMap<>(5*nKmers/3);
        for ( Tuple2<Kmer,Set<Long>> tupl : readPullingKmers )
            readPullingMap.put(tupl._1,new ArrayList<>(tupl._2));

        return readPullingMap;
    }

    // writes a stream of reads to a FASTQ file
    private static void writeFastq( Long breakpointId, Iterator<GATKRead> readItr, String outputDir )
    {
        File fastqName = new File(outputDir, "assembly" + breakpointId + ".fastq");
        try
        {
            BufferedWriter writer = new BufferedWriter(new FileWriter(fastqName));
            while ( readItr.hasNext() )
            {
                GATKRead read = readItr.next();
                writer.write("@");
                writer.write(read.getName());
                writer.write(read.isFirstOfPair()?":1":":2");
                writer.newLine();
                writer.write(read.getBasesString());
                writer.newLine();
                writer.write('+');
                writer.newLine();
                for ( byte b : read.getBaseQualities() )
                    writer.write(33+b);
                writer.newLine();
            }
            writer.close();
        }
        catch ( IOException e )
        {
            throw new GATKException("Can't write "+fastqName+'.',e);
        }
    }

    private static final class EventLocus implements Comparable<EventLocus>
    {
        private final int locusStart;
        private final int locusEnd;
        private long id = 0; // used as a uniquifier

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
    }

    // Class to detect a split read:  that's a read where one part of the read maps well to reference, but
    // another portion maps elsewhere or not at all.
    // It returns the approximate locus of the breakpoint if it can find one.
    // Class assumes that all reads presented to it are mapped.
    private static final class SplitReadDetector
    {
        private static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
        private static final int MIN_INDEL_LEN = 25; // minimum length of an interesting indel
        private static final int SOFT_CLIP_LOCUS_WIDTH = 4; // uncertainty in event locus for soft clip
        private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window

        public EventLocus getLocusIfReadIsSplit( GATKRead read )
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

            int locus = read.getStart();
            for ( CigarElement ele : cigarElements )
            {
                CigarOperator op = ele.getOperator();
                if ( op == CigarOperator.DELETION || op == CigarOperator.INSERTION )
                    if ( ele.getLength() >= MIN_INDEL_LEN )
                        return new EventLocus(locus-SOFT_CLIP_LOCUS_WIDTH/2, SOFT_CLIP_LOCUS_WIDTH);
                if ( op.consumesReferenceBases() )
                    locus += ele.getLength();
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
    }

    // Class to detect reads that are part of a discordant pair.
    // Concordant pairs are mapped to opposite strands of a reference sequence at a distance that is consistant with the
    // library fragment length.
    // Discordant pairs aren't concordant, except that we're also ignoring reads with mates mapped to other contigs.
    // Class assumes that all reads presented to it are mapped.
    private static final class DiscordantPairDetector
    {
        private static final int FUNKY_PAIR_LOCUS_WIDTH = 100; // uncertainty in event locus for discordant pair
        //TODO: Above should depend on the fragment size distribution

        public EventLocus getLocusIfDiscordant( GATKRead read )
        {
            if ( read.mateIsUnmapped() )
            {
                int locus = read.isReverseStrand() ? read.getStart() - FUNKY_PAIR_LOCUS_WIDTH : read.getEnd();
                return new EventLocus(locus,FUNKY_PAIR_LOCUS_WIDTH);
            }
            else if ( read.getContig() == read.getMateContig() )
            {
                if ( read.isReverseStrand() == read.mateIsReverseStrand() || !read.isProperlyPaired() )
                {
                    int locus = read.getFragmentLength() < 0 ? read.getStart() - FUNKY_PAIR_LOCUS_WIDTH : read.getEnd();
                    return new EventLocus(locus,FUNKY_PAIR_LOCUS_WIDTH);
                }
            }

            return null;
        }
    }

    // A class that acts as a filter for a stream of reads.  It passes only those that are a part of a putative
    // breakpoint cluster.
    private static final class BreakpointClusterer implements Iterator<GATKRead>, Iterable<GATKRead>
    {
        private SplitReadDetector splitReadDetector;
        private DiscordantPairDetector discordantPairDetector;
        private final Iterator<GATKRead> inputIterator;
        private String currentContig = null;
        private long eventIdCounter = 0;
        private final SortedMap<EventLocus,GATKRead> locMap = new TreeMap<>();
        private Iterator<GATKRead> outputIterator = Collections.emptyIterator();

        private static final int MIN_MATCH_LEN = 45; // minimum length of matched portion of an interesting alignment
        private static final int MIN_MAPQ = 20; // minimum mapping quality of an interesting alignment
        private static final int MAX_LOCUS_DIST = 300; // stale evidence distance, should be about the fragment length
        private static final int MIN_EVIDENCE = 5; // minimum evidence count in a cluster

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

        @Override
        public Iterator<GATKRead> iterator() { return this; }

        private Iterator<GATKRead> processRead( final GATKRead read )
        {
            if ( !read.isUnmapped() &&
                    read.getMappingQuality() >= MIN_MAPQ &&
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
                    EventLocus eventLocus = splitReadDetector.getLocusIfReadIsSplit(read);
                    if ( eventLocus != null )
                        return cluster(eventLocus,read);
                    eventLocus = discordantPairDetector.getLocusIfDiscordant(read);
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
    }

    // Class that acts as a mapper from a stream of reads to a stream of <Kmer,BreakpointId> pairs.
    // Sitting atop the BreakpointClusterer, it groups funky reads that are mapped near each other and assigns that
    // group a breakpoint id, and transforms each funky read into kmers which are emitted along with the breakpoint id.
    private static final class ClusterKmerizer implements Iterator<Tuple2<Kmer,Long>>, Iterable<Tuple2<Kmer,Long>>
    {
        private final Iterator<GATKRead> inputIterator;
        private final Set<Kmer> killSet;
        private Iterator<Tuple2<Kmer,Long>> outputIterator = Collections.emptyIterator();
        private String currentContig = null;
        private int lastLocus = 0;
        private Long breakpointId;

        private static final int CLUSTER_BREAK_DISTANCE = 500;
        //TODO: Should depend on fragment size

        public ClusterKmerizer( Integer partIdx, Iterator<GATKRead> someInputIterator, Set<Kmer> someKillSet )
        {
            inputIterator = someInputIterator;
            killSet = someKillSet;
            breakpointId = partIdx.longValue() << 32;
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
        public Tuple2<Kmer,Long> next()
        {
            if ( !hasNext() )
                throw new NoSuchElementException("ClusterKmerizer is exhausted.");
            return outputIterator.next();
        }

        @Override
        public Iterator<Tuple2<Kmer,Long>> iterator() { return this; }

        private Iterator<Tuple2<Kmer,Long>> processRead( GATKRead read )
        {
            if ( read.getContig() != currentContig || read.getStart() > lastLocus+CLUSTER_BREAK_DISTANCE )
                breakpointId += 1;
            currentContig = read.getContig();
            lastLocus = read.getStart();

            List<Tuple2<Kmer,Long>> kmerList = new ArrayList<>(read.getLength()-KLEN+1);
            BytesKmerizer kmerizer = new BytesKmerizer(read.getBases(),KLEN);
            while ( kmerizer.hasNext() )
            {
                Kmer kmer = kmerizer.next().canonical(KLEN);
                if ( !killSet.contains(kmer) )
                    kmerList.add(new Tuple2<>(kmer,breakpointId));
            }
            return kmerList.iterator();
        }
    }

    // Class that acts as a mapper from a stream of reads to a stream of <breakpointId,read> pairs.
    // It knows which breakpoint(s) a read belongs to (if any) by kmerizing the read, and looking up each Kmer in a
    // map of Kmers onto breakpoints.
    private static final class ReadPuller
            implements Iterator<Tuple2<Long,GATKRead>>, Iterable<Tuple2<Long,GATKRead>>
    {
        private final Iterator<GATKRead> inputIterator;
        private final Map<Kmer,List<Long>> readPullingMap;
        private Iterator<Tuple2<Long,GATKRead>> outputIterator = Collections.emptyIterator();

        public ReadPuller( Iterator<GATKRead> someInputIterator, Map<Kmer,List<Long>> someReadPullingMap )
        {
            inputIterator = someInputIterator;
            readPullingMap = someReadPullingMap;
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
        public Tuple2<Long,GATKRead> next()
        {
            if ( !hasNext() )
                throw new NoSuchElementException("ClusterKmerizer is exhausted.");
            return outputIterator.next();
        }

        @Override
        public Iterator<Tuple2<Long,GATKRead>> iterator() { return this; }

        private Iterator<Tuple2<Long,GATKRead>> processRead( GATKRead read )
        {
            Set<Long> breakpointSet = new HashSet<>();
            BytesKmerizer kmerizer = new BytesKmerizer(read.getBases(),KLEN);
            while ( kmerizer.hasNext() )
            {
                Kmer kmer = kmerizer.next().canonical(KLEN);
                List<Long> breakpoints = readPullingMap.get(kmer);
                if ( breakpoints != null )
                    breakpointSet.addAll(breakpoints);
            }

            if ( breakpointSet.isEmpty() )
                return Collections.emptyIterator();

            List<Tuple2<Long,GATKRead>> readList = new ArrayList<>(breakpointSet.size());
            for ( Long breakpointId : breakpointSet )
                readList.add(new Tuple2<>(breakpointId,read));

            return readList.iterator();
        }
    }
/*
    // dead code to reorder a nearly coordinate-ordered stream of reads into their original coordinate order
    // i'd like to keep this in here for one PR, just so it's part of the "permanent record", and then i'll delete it

    private static final class WindowSorter implements Iterator<GATKRead>, Iterable<GATKRead>
    {
        private String currentContig;
        private int contigNo;
        private final SortedSet<GATKKey> recordSet;
        private final Iterator<GATKRead> inputIterator;
        private final int windowSize;

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

        private static final class GATKKey implements Comparable<GATKKey>
        {
            private final int relativeContigNo;
            private final int locus;
            private final boolean isFirst;
            private final GATKRead read;

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
        }
    }
*/
}
