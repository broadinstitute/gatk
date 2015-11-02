package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.MapPartitioner;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadConstants;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * SparkTool to produce small FASTQs of reads sharing kmers with putative SV breakpoints for local assembly.
 */
@CommandLineProgramProperties(summary="Produce FASTQs of breakpoint candidates",
        oneLineSummary="Produce small FASTQs of reads sharing kmers with putative SV breakpoints for local assembly",
        programGroup = SparkProgramGroup.class)
public class FindSVBreakpointsSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "dir for output FASTQs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputDir;

    /**
     * This is a path to a file of kmers that appear too frequently in the reference to be usable as probes to localize
     * reads.  We don't calculate it here, because it depends only on the reference.
     * The program FindBadGenomicKmers can produce such a list for you.
     */
    @Argument(doc = "file containing ubiquitous kmer list", shortName = "KS", fullName = "kmersToIgnore",
            optional = false)
    private String kmersToIgnoreFilename;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader header = getHeaderForReads();
        if ( header.getSortOrder() != SAMFileHeader.SortOrder.coordinate ) {
            throw new GATKException("The BAM must be coordinate sorted.");
        }

        // Process the input data to produce a map of SVKmer -> List of breakpoint IDs.
        // (A breakpoint ID is just an arbitrary identifier for the putative breakpoint.)
        // See the doc for findClusters to see a description of how this big Spark pipeline works.
        final Broadcast<Map<SVKmer, List<Long>>> broadcastMap = ctx.broadcast(findClusters(ctx));

        // Process the input again, this time pulling all reads that contain kmers associated with a given breakpoint,
        // and writing those reads into a separate FASTQ for each breakpoint.
        getUnfilteredReads()
            .filter(read ->
                    !read.isSecondaryAlignment() && !read.isSupplementaryAlignment() &&
                    !read.isDuplicate() && !read.failsVendorQualityCheck())
            .mapPartitionsToPair(readItr ->
                    new MapPartitioner<>(readItr,new ReadsForBreakpointFinder(broadcastMap.value())))
            .groupByKey()
            .foreach(tupl -> writeFastq(tupl._1, tupl._2.iterator(), outputDir));
    }

    /**
     * Read a file of kmers to ignore.
     * Each line must be exactly SVConstants.KMER_SIZE characters long, and must match [ACGT]*.
     */
    @VisibleForTesting static Set<SVKmer> readKmersToIgnore( final File kmersToIgnoreFile ) {
        // compute size of HashMap to provide a final load factor of 7/10 (a bit lower than default max which is 3/4)
        final int nKmers = 10*(int)(kmersToIgnoreFile.length()/(SVConstants.KMER_SIZE+1))/7+1;
        final Set<SVKmer> kmersToIgnore = new HashSet<>(nKmers);

        try ( final BufferedReader rdr = new BufferedReader(new FileReader(kmersToIgnoreFile)) ) {
            String line;
            while ( (line = rdr.readLine()) != null ) {
                if ( line.length() != SVConstants.KMER_SIZE ) {
                    throw new GATKException("SVKmer kill set contains a line of length " + line.length() +
                            " but we were expecting K=" + SVConstants.KMER_SIZE);
                }

                final SVKmerizer kmerizer = new SVKmerizer(line, SVConstants.KMER_SIZE);
                if ( !kmerizer.hasNext() ) {
                    throw new GATKException("Unable to kmerize the kmer kill set string '" + line + "'.");
                }

                kmersToIgnore.add(kmerizer.next());
            }
        }
        catch ( final IOException e ) {
            throw new GATKException("Unable to read kmer kill set.", e);
        }

        return kmersToIgnore;
    }

    /**
     * Produce a map of SVKmer -> List of unique breakpoint IDs.
     * This is done by looking for reads that appear to be split, or for reads that are part of discordant pairs.
     * When a sufficient number of these reads are mapped onto reference near one another, those reads are taken
     * together as evidence of a putative breakpoint.
     * Each such breakpoint is numbered, and the reads that supported that breakpoint are kmerized.
     * Kmers that appear on a kill-list are ignored.
     * This produces a map of Breakpoint ID -> List of Kmers, which is inverted to produce the final output.
     */
    private Map<SVKmer, List<Long>> findClusters( final JavaSparkContext ctx ) {
        // read the ignored kmer list
        final Broadcast<Set<SVKmer>> kmersToIgnore = ctx.broadcast(readKmersToIgnore(new File(kmersToIgnoreFilename)));

        // produce a set of Breakpoint IDs for each of an interesting set of Kmers by...
        // 1) identifying putative breakpoints supported by a cluster of funky reads
        // 2) kmerizing all the funky reads in the cluster
        // 3) grouping by SVKmer
        final List<Tuple2<SVKmer, List<Long>>> readPullingKmers =
            getUnfilteredReads()
                .mapPartitionsWithIndex((partitionIdx, readItr) ->
                    new MapPartitioner<>(new MapPartitioner<>(readItr, new BreakpointClusterer()).iterator(),
                            new ClusterKmerizer(partitionIdx,kmersToIgnore.value())).iterator(),
                    false )
                .mapToPair(tupl -> tupl)
                .combineByKey(breakId -> { List<Long> ids = new ArrayList<>(1); ids.add(breakId); return ids; },
                        (ids, breakId) -> { if ( !ids.contains(breakId) ) ids.add(breakId); return ids; },
                        (ids1, ids2) -> {
                            for ( Long breakId : ids2 ) if ( !ids1.contains(breakId) ) ids1.add(breakId); return ids1;
                        })
                .collect();

        // turn the list of interesting Kmers and the set of putative breakpoints to which they belong into a map

        // compute size of HashMap to provide a final load factor of 7/10 (a bit lower than default max which is 3/4)
        final int hashMapSize = 10*readPullingKmers.size()/7 + 1;
        final Map<SVKmer, List<Long>> readPullingMap = new HashMap<>(hashMapSize);
        for ( final Tuple2<SVKmer, List<Long>> tupl : readPullingKmers ) {
            readPullingMap.put(tupl._1, tupl._2);
        }

        return readPullingMap;
    }

    /** writes a stream of reads to a FASTQ file */
    @VisibleForTesting static File writeFastq( final Long breakpointId, final Iterator<GATKRead> readItr, final String outputDir ) {
        final File fastqName = new File(outputDir, "assembly" + breakpointId + ".fastq");
        try ( final FastqWriter writer = new FastqWriterFactory().newWriter(fastqName) ) {
            while ( readItr.hasNext() ) {
                final GATKRead read = readItr.next();
                String readName = read.getName();
                if ( read.isPaired() ) {
                    readName += read.isSecondOfPair() ? "/2" : "/1";
                }
                final String quals = ReadUtils.getBaseQualityString(read);
                writer.write(new FastqRecord(readName, read.getBasesString(), "", quals));
            }
        }
        return fastqName;
    }

    @VisibleForTesting static final class EventLocus implements Comparable<EventLocus> {
        private final int locusStart;
        private final int locusEnd;
        private long id = 0; // used as a uniquifier for lack of a multi-map

        public EventLocus( final int locusStart, final int intervalWidth, final long id ) {
            this.locusStart = locusStart;
            this.locusEnd = locusStart + intervalWidth;
            this.id = id;
        }

        public int getLocusStart() { return locusStart; }
        public int getLocusEnd() { return locusEnd; }

        @Override
        public int compareTo( final EventLocus that ) {
            int result = Integer.compare(locusStart, that.locusStart);
            if ( result == 0 ) result = Integer.compare(locusEnd, that.locusEnd);
            if ( result == 0 ) result = Long.compare(id, that.id);
            return result;
        }
    }

    /**
     * Class to detect a split read:  that's a read where one part of the read maps well to reference, but
     * another portion maps elsewhere or not at all.
     * It returns the approximate locus of the breakpoint if it can find one.
     * Class assumes that all reads presented to it are mapped.
     */
    @VisibleForTesting static final class SplitReadDetector {
        private static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
        private static final int MIN_INDEL_LEN = 25; // minimum length of an interesting indel
        private static final int SOFT_CLIP_LOCUS_WIDTH = 4; // uncertainty in event locus for soft clip
        private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window

        public EventLocus getLocusIfReadIsSplit( final GATKRead read, final long readId ) {
            final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
            if ( hasInitialSoftClip(cigarElements,read) ) {
                return new EventLocus(read.getStart() - SOFT_CLIP_LOCUS_WIDTH / 2, SOFT_CLIP_LOCUS_WIDTH, readId);
            }

            if ( hasFinalSoftClip(cigarElements,read) ) {
                return new EventLocus(read.getEnd() - SOFT_CLIP_LOCUS_WIDTH / 2, SOFT_CLIP_LOCUS_WIDTH, readId);
            }

            // or if the read has a big indel, return that as the locus
            int locus = read.getStart();
            for ( final CigarElement ele : cigarElements ) {
                final CigarOperator op = ele.getOperator();
                if ( op == CigarOperator.DELETION || op == CigarOperator.INSERTION ) {
                    if ( ele.getLength() >= MIN_INDEL_LEN ) {
                        return new EventLocus(locus - SOFT_CLIP_LOCUS_WIDTH / 2, SOFT_CLIP_LOCUS_WIDTH, readId);
                    }
                }
                if ( op.consumesReferenceBases() )
                    locus += ele.getLength();
            }
            return null;
        }

        private static boolean hasInitialSoftClip( final List<CigarElement> cigarElements, final GATKRead read ) {
            final ListIterator<CigarElement> itr = cigarElements.listIterator();
            if ( !itr.hasNext() ) return false;

            CigarElement firstEle = itr.next();
            if ( firstEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasNext() ) {
                firstEle = itr.next();
            }
            return firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                    firstEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                    isHighQualityRegion(read.getBaseQualities(), 0);
        }

        private static boolean hasFinalSoftClip( final List<CigarElement> cigarElements, final GATKRead read ) {
            final ListIterator<CigarElement> itr = cigarElements.listIterator(cigarElements.size());
            if ( !itr.hasPrevious() ) return false;

            CigarElement lastEle = itr.previous();
            if ( lastEle.getOperator() == CigarOperator.HARD_CLIP && itr.hasPrevious() ) {
                lastEle = itr.previous();
            }
            return lastEle.getOperator() == CigarOperator.SOFT_CLIP &&
                    lastEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                    isHighQualityRegion(read.getBaseQualities(), read.getLength() - lastEle.getLength());
        }

        private static boolean isHighQualityRegion( final byte[] quals, int idx ) {
            for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx ) {
                if ( quals[idx] < MIN_QUALITY ) return false;
            }
            return true;
        }
    }

    /**
     * Class to detect reads that are part of a discordant pair.
     * Concordant pairs are mapped to opposite strands of a reference sequence at a distance that is consistant with
     * the library fragment length.
     * Discordant pairs aren't concordant, except that we're also ignoring reads with mates mapped to other contigs.
     * Class assumes that all reads presented to it are mapped.
     */
    @VisibleForTesting static final class DiscordantPairDetector {
        private static final int FUNKY_PAIR_LOCUS_WIDTH = 100; // uncertainty in event locus for discordant pair
        //TODO: Above should depend on the fragment size distribution

        public EventLocus getLocusIfDiscordant( final GATKRead read, final long readId ) {
            if ( read.mateIsUnmapped() ) {
                final int locus = read.isReverseStrand() ? read.getStart() - FUNKY_PAIR_LOCUS_WIDTH : read.getEnd();
                return new EventLocus(locus, FUNKY_PAIR_LOCUS_WIDTH, readId);
            } else if ( read.getContig().equals(read.getMateContig()) ) {
                if ( read.isReverseStrand() == read.mateIsReverseStrand() || !read.isProperlyPaired() ) {
                    final int locus = read.getFragmentLength() < 0 ? read.getStart() - FUNKY_PAIR_LOCUS_WIDTH : read.getEnd();
                    return new EventLocus(locus, FUNKY_PAIR_LOCUS_WIDTH, readId);
                }
            }

            return null;
        }
    }

    /**
     * A class that acts as a filter for a stream of reads.  It passes only those that are a part of a putative
     * breakpoint cluster.
     */
    @VisibleForTesting static final class BreakpointClusterer implements Function<GATKRead, Iterator<GATKRead>> {
        private final SplitReadDetector splitReadDetector;
        private final DiscordantPairDetector discordantPairDetector;
        private String currentContig = null;
        private long readId = 0;
        private final SortedMap<EventLocus, GATKRead> locMap = new TreeMap<>();

        private static final int MIN_MATCH_LEN = 45; // minimum length of matched portion of an interesting alignment
        @VisibleForTesting static final int MIN_MAPQ = 20; // minimum mapping quality of an interesting alignment
        private static final int MAX_LOCUS_DIST = 300; // stale evidence distance, should be about the fragment length
        @VisibleForTesting static final int MIN_EVIDENCE = 5; // minimum evidence count in a cluster

        public BreakpointClusterer() {
            this.splitReadDetector = new SplitReadDetector();
            this.discordantPairDetector = new DiscordantPairDetector();
        }

        public Iterator<GATKRead> apply( final GATKRead read ) {
            readId += 1;
            if ( !read.isSecondaryAlignment() && !read.isSupplementaryAlignment() &&
                 !read.isDuplicate() && !read.failsVendorQualityCheck() && !read.isUnmapped() &&
                 read.getMappingQuality() >= MIN_MAPQ &&
                 read.getStart() != ReadConstants.UNSET_POSITION ) {
                final String readContig = read.getContig();
                if ( currentContig == null || !currentContig.equals(readContig) ) {
                    currentContig = readContig;
                    locMap.clear();
                }

                final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
                int matchLen = 0;
                for ( final CigarElement ele : cigarElements ) {
                    if ( ele.getOperator() == CigarOperator.MATCH_OR_MISMATCH ) {
                        matchLen += ele.getLength();
                    }
                }
                if ( matchLen >= MIN_MATCH_LEN ) {
                    EventLocus eventLocus = splitReadDetector.getLocusIfReadIsSplit(read, readId);
                    if ( eventLocus != null ) return cluster(eventLocus, read);
                    eventLocus = discordantPairDetector.getLocusIfDiscordant(read, readId);
                    if ( eventLocus != null ) return cluster(eventLocus, read);
                }
            }
            return Collections.emptyIterator();
        }

        private Iterator<GATKRead> cluster( final EventLocus newEventLocus, final GATKRead read ) {
            final int locusStart = newEventLocus.getLocusStart();
            final int locusEnd = newEventLocus.getLocusEnd();

            // clean out old stuff that can't possibly be interesting anymore
            removeStaleMapEntries(locusStart, locMap.entrySet().iterator());

            locMap.put(newEventLocus, read);

            // see if there are a sufficient number of overlapping putative events to call a cluster
            final List<Map.Entry<EventLocus, GATKRead>> overlappers = new LinkedList<>();
            for ( final Map.Entry<EventLocus, GATKRead> entry : locMap.entrySet() ) {
                final EventLocus eventLocus = entry.getKey();
                if ( eventLocus.getLocusStart() >= locusEnd ) break;
                if ( eventLocus.getLocusEnd() > locusStart ) {
                    overlappers.add(entry);
                }
            }
            if ( overlappers.size() >= MIN_EVIDENCE ) {
                final List<GATKRead> readsList = new LinkedList<>();
                overlappers.stream()
                        .filter(entry -> entry.getValue() != null)
                        .forEach(entry -> { readsList.add(entry.getValue()); entry.setValue(null); });
                return readsList.iterator();
            }
            return Collections.emptyIterator();
        }

        private static void removeStaleMapEntries( final int locusStart,
                                                   final Iterator<Map.Entry<EventLocus, GATKRead>> itr ) {
            final int staleEnd = locusStart - MAX_LOCUS_DIST;
            while ( itr.hasNext() ) {
                final EventLocus eventLocus = itr.next().getKey();
                if ( eventLocus.getLocusStart() >= locusStart ) break;
                if ( eventLocus.getLocusEnd() <= staleEnd ) {
                    itr.remove();
                }
            }
        }
    }

    /**
     * Class that acts as a mapper from a stream of reads to a stream of <SVKmer,BreakpointId> pairs.
     * Sitting atop the BreakpointClusterer, it groups funky reads that are mapped near each other and assigns that
     * group a breakpoint id, and transforms each funky read into kmers which are emitted along with the breakpoint id.
     */
    @VisibleForTesting static final class ClusterKmerizer implements Function<GATKRead, Iterator<Tuple2<SVKmer, Long>>> {
        private final Set<SVKmer> kmersToIgnore;
        private String currentContig = null;
        private int lastLocus = 0;
        private Long breakpointId;

        private static final int CLUSTER_BREAK_DISTANCE = 500;
        //TODO: Should depend on fragment size

        public ClusterKmerizer( final Integer partitionIdx, final Set<SVKmer> kmersToIgnore ) {
            this.kmersToIgnore = kmersToIgnore;
            // to make breakpointIds unique across all partitions of the input data, we'll
            // start each ClusterKmerizer with a breakpointId of partitionIdx * 2^32.
            this.breakpointId = partitionIdx.longValue() << Integer.SIZE;
        }

        public Iterator<Tuple2<SVKmer, Long>> apply( final GATKRead read ) {
            if ( !read.getContig().equals(currentContig) || read.getStart() > lastLocus+CLUSTER_BREAK_DISTANCE ) {
                breakpointId += 1;
            }
            currentContig = read.getContig();
            lastLocus = read.getStart();

            final int nKmers = read.getLength() - SVConstants.KMER_SIZE + 1;
            return SVKmerizer.stream(read.getBases(), SVConstants.KMER_SIZE)
                .map(kmer -> kmer.canonical(SVConstants.KMER_SIZE))
                .filter(kmer -> !kmersToIgnore.contains(kmer))
                .map(kmer -> new Tuple2<>(kmer,breakpointId))
                .collect(Collectors.toCollection(() -> new ArrayList<>(nKmers)))
                .iterator();
        }
    }

    /**
     * Class that acts as a mapper from a stream of reads to a stream of <breakpointId,read> pairs.
     * It knows which breakpoint(s) a read belongs to (if any) by kmerizing the read, and looking up each SVKmer in
     * a map of Kmers onto breakpoints.
     */
    @VisibleForTesting static final class ReadsForBreakpointFinder
            implements Function<GATKRead, Iterator<Tuple2<Long, GATKRead>>> {
        private final Map<SVKmer, List<Long>> readPullingMap;

        public ReadsForBreakpointFinder( final Map<SVKmer, List<Long>> readPullingMap ) {
            this.readPullingMap = readPullingMap;
        }

        public Iterator<Tuple2<Long, GATKRead>> apply( final GATKRead read ) {
            final HashSet<Long> breakpoints = SVKmerizer.stream(read.getBases(), SVConstants.KMER_SIZE)
                .map(kmer -> kmer.canonical(SVConstants.KMER_SIZE))
                .flatMap(kmer -> readPullingMap.getOrDefault(kmer, Collections.emptyList()).stream())
                .collect(Collectors.toCollection(HashSet::new));
            final int nBreakpoints = breakpoints.size();
            return breakpoints
                .stream()
                .map(breakpointId -> new Tuple2<>(breakpointId, read))
                .collect(Collectors.toCollection(() -> new ArrayList<>(nBreakpoints)))
                .iterator();
        }
    }
}
