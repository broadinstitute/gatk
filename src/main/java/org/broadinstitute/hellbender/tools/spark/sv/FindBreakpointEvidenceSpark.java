package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import org.apache.commons.collections4.iterators.SingletonIterator;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tool to describe reads that support a hypothesis of a genomic breakpoint.
 */
@CommandLineProgramProperties(summary="Find reads that evidence breakpoints.",
        oneLineSummary="Dump FASTQs for local assembly of putative genomic breakpoints.",
        programGroup = StructuralVariationSparkProgramGroup.class)
public final class FindBreakpointEvidenceSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    //--------- parameters ----------

    // no-arg constructor for Params object establishes default values
    @VisibleForTesting static final Params defaultParams = new Params();

    @Argument(doc = "Kmer size.", fullName = "kSize")
    private int kSize = defaultParams.kSize;

    @Argument(doc = "Minimum kmer entropy", fullName = "kmerEntropy")
    private double minEntropy = defaultParams.minEntropy;

    @Argument(doc = "The minimum mapping quality for reads used to gather evidence of breakpoints.",
            fullName = "minEvidenceMapQ", optional = true)
    private int minEvidenceMapQ = defaultParams.minEvidenceMapQ;

    @Argument(doc = "The minimum length of the matched portion of an interesting alignment.  "+
            "Reads that don't match at least this many reference bases won't be used in gathering evidence.",
            fullName = "minEvidenceMatchLength", optional = true)
    private int minEvidenceMatchLength = defaultParams.minEvidenceMatchLength;

    @Argument(doc = "Intervals with more than this much coverage are filtered out, because the reads mapped to "+
            "that interval are clearly not exclusively local to the interval.", fullName = "maxIntervalCoverage")
    private int maxIntervalCoverage = defaultParams.maxIntervalCoverage;

    @Argument(doc = "Minimum number of reads in cluster to declare an interval of interest.",
            fullName = "minEvidenceCount")
    private int minEvidenceCount = defaultParams.minEvidenceCount;

    @Argument(doc = "Minimum number of localizing kmers in a valid interval.", fullName="minKmersPerInterval")
    private int minKmersPerInterval = defaultParams.minKmersPerInterval;

    @Argument(doc = "KmerCleaner maximum number of intervals for a localizing kmer.", fullName = "cleanerMaxIntervals")
    private int cleanerMaxIntervals = defaultParams.cleanerMaxIntervals;

    @Argument(doc = "KmerCleaner minimum kmer count.", fullName = "cleanerMinKmerCount")
    private int cleanerMinKmerCount = defaultParams.cleanerMinKmerCount;

    @Argument(doc = "KmerCleaner maximum kmer count.", fullName = "cleanerMaxKmerCount")
    private int cleanerMaxKmerCount = defaultParams.cleanerMaxKmerCount;

    @Argument(doc = "KmerCleaner unique error-free kmers per partition", fullName = "cleanerKmersPerPartitionGuess")
    private int cleanerKmersPerPartitionGuess = defaultParams.cleanerKmersPerPartitionGuess;

    @Argument(doc = "Maximum number of templates containing an assembly kmer.", fullName = "maxQNamesPerKmer")
    private int maxQNamesPerKmer = defaultParams.maxQNamesPerKmer;

    @Argument(doc = "Guess at number of clean kmers per assembly partition.", fullName = "assemblyKmerMapSize")
    private int assemblyKmerMapSize = defaultParams.assemblyKmerMapSize;

    @Argument(doc = "Guess at the ratio of reads in the final assembly to the number reads mapped to the interval.",
            fullName = "assemblyToMappedSizeRatioGuess")
    private int assemblyToMappedSizeRatioGuess = defaultParams.assemblyToMappedSizeRatioGuess;

    @Argument(doc = "Maximum FASTQ file size.", fullName = "maxFASTQSize")
    private int maxFASTQSize = defaultParams.maxFASTQSize;

    @Argument(doc = "Exclusion interval padding.", fullName = "exclusionIntervalPadding")
    private int exclusionIntervalPadding = defaultParams.exclusionIntervalPadding;

    @Argument(doc = "Include read mapping location in FASTQ files.", fullName = "includeMappingLocation")
    private boolean includeMappingLocation = true;

    @Argument(doc = "Include read mapping location in FASTQ files.", fullName = "intervalOnlyAssembly")
    private boolean intervalOnlyAssembly = false;

    // --------- locations ----------

    @Argument(doc = "directory for fastq output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

    @Argument(doc = "file for read metadata", fullName = "readMetadata", optional = true)
    private String metadataFile;

    @Argument(doc = "directory for evidence output", fullName = "breakpointEvidenceDir", optional = true)
    private String evidenceDir;

    @Argument(doc = "file for breakpoint intervals output", fullName = "breakpointIntervals", optional = true)
    private String intervalFile;

    @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsMapped", optional = true)
    private String qNamesMappedFile;

    @Argument(doc = "file for kmer intervals output", fullName = "kmerIntervals", optional = true)
    private String kmerFile;

    @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsForAssembly", optional = true)
    private String qNamesAssemblyFile;

    /**
     * This is a file that calls out the coordinates of intervals in the reference assembly to exclude from
     * consideration when calling putative breakpoints.
     * Each line is a tab-delimited interval with 1-based inclusive coordinates like this:
     *  chr1	124535434	142535434
     */
    @Argument(doc = "file of reference intervals to exclude", fullName = "exclusionIntervals", optional = true)
    private String exclusionIntervalsFile;

    /**
     * This is a path to a file of kmers that appear too frequently in the reference to be usable as probes to localize
     * reads.  We don't calculate it here, because it depends only on the reference.
     * The program FindBadGenomicKmersSpark can produce such a list for you.
     */
    @Argument(doc = "file containing ubiquitous kmer list. see FindBadGenomicKmersSpark to generate it.",
            fullName = "kmersToIgnore")
    private String kmersToIgnoreFile;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        final SAMFileHeader header = getHeaderForReads();
        if ( header.getSortOrder() != SAMFileHeader.SortOrder.coordinate ) {
            throw new GATKException("The reads must be coordinate sorted.");
        }

        final Locations locations =
                new Locations(metadataFile, evidenceDir, intervalFile, qNamesMappedFile,
                                kmerFile, qNamesAssemblyFile, exclusionIntervalsFile);
        final Params params =
                new Params(kSize, minEntropy, minEvidenceMapQ, minEvidenceMatchLength, maxIntervalCoverage,
                            minEvidenceCount, minKmersPerInterval, cleanerMaxIntervals, cleanerMinKmerCount,
                            cleanerMaxKmerCount, cleanerKmersPerPartitionGuess, maxQNamesPerKmer, assemblyKmerMapSize,
                            assemblyToMappedSizeRatioGuess, maxFASTQSize, exclusionIntervalPadding);

        final PipelineOptions pipelineOptions = getAuthenticatedGCSOptions();
        final JavaRDD<GATKRead> unfilteredReads = getUnfilteredReads();
        final JavaRDD<GATKRead> allPrimaryLines =
                unfilteredReads.filter(read -> !read.isSecondaryAlignment() && !read.isSupplementaryAlignment());

        // develop evidence, intervals, and, finally, a set of template names for each interval
        final Tuple2<List<SVInterval>, HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> intervalsAndQNameMap =
                getMappedQNamesSet(params, ctx, header, unfilteredReads, locations, pipelineOptions);
        final List<SVInterval> intervals = intervalsAndQNameMap._1;
        if ( intervals.isEmpty() ) return;

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap = intervalsAndQNameMap._2;

        // supplement the template names with other reads that share kmers
        final Map<Integer, String> intervalDispositions;
        if ( intervalOnlyAssembly ) {
            intervalDispositions = new HashMap<>();
        } else {
            intervalDispositions = addAssemblyQNames(params, ctx, kmersToIgnoreFile, qNamesMultiMap, intervals.size(),
                    allPrimaryLines, locations, pipelineOptions);
        }

        // write a FASTQ file for each interval
        final String outDir = outputDir;
        final int maxFastqSize = maxFASTQSize;
        final boolean includeMapLoc = includeMappingLocation;
        intervalDispositions.putAll(
                generateFastqs(ctx, qNamesMultiMap, allPrimaryLines, intervals.size(), includeMapLoc,
                                intervalAndFastqBytes -> writeFastq(intervalAndFastqBytes, outDir, maxFastqSize)));

        // record the intervals
        if ( locations.intervalFile != null ) {
            writeIntervalFile(locations.intervalFile, pipelineOptions, header, intervals, intervalDispositions);
        }

        log("Wrote FASTQs for assembly.");
    }

    /** write a file describing each interval */
    private static void writeIntervalFile( final String intervalFile,
                                           final PipelineOptions pipelineOptions,
                                           final SAMFileHeader header,
                                           final List<SVInterval> intervals,
                                           final Map<Integer, String> intervalDispositions ) {

        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(intervalFile, pipelineOptions)))) {
            final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
            final int nIntervals = intervals.size();
            for ( int intervalId = 0; intervalId != nIntervals; ++intervalId ) {
                final SVInterval interval = intervals.get(intervalId);
                final String seqName = contigs.get(interval.getContig()).getSequenceName();
                String disposition = intervalDispositions.get(intervalId);
                if ( disposition == null ) disposition = "unknown";
                writer.write(intervalId + "\t" +
                        seqName + ":" + interval.getStart() + "-" + interval.getEnd() + "\t" +
                        disposition + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write intervals file " + intervalFile, ioe);
        }
    }

    /**
     * Find the breakpoint evidence,
     * cluster the evidence into intervals,
     * find the template names mapped into each interval,
     * kmerize each of these templates,
     * clean up by removing some intervals that are bogus as evidenced by ubiquitous kmers,
     * and return a set of template names and the intervals to which they belong.
     */
    private Tuple2<List<SVInterval>, HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> getMappedQNamesSet(
            final Params params,
            final JavaSparkContext ctx,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final JavaRDD<GATKRead> mappedReads =
                unfilteredReads.filter(read ->
                        !read.isDuplicate() && !read.failsVendorQualityCheck() && !read.isUnmapped());
        final ReadMetadata readMetadata = new ReadMetadata(header, mappedReads);
        if ( locations.metadataFile != null ) writeMetadata(readMetadata, locations.metadataFile, pipelineOptions);
        log("Metadata retrieved.");

        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);
        List<SVInterval> intervals = getIntervals(params, broadcastMetadata, header, mappedReads, locations);

        final int nIntervals = intervals.size();
        log("Discovered " + nIntervals + " intervals.");

        if ( nIntervals == 0 ) return null;

        if ( locations.exclusionIntervalsFile != null ) {
            intervals = removeIntervalsNearGapsAndLog(intervals, params.exclusionIntervalPadding, readMetadata,
                    locations.exclusionIntervalsFile, pipelineOptions);
        }
        intervals = removeHighCoverageIntervalsAndLog(params, ctx, broadcastMetadata, intervals, mappedReads);

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap =
                getQNames(params, ctx, broadcastMetadata, intervals, mappedReads);
        broadcastMetadata.destroy();

        if ( locations.qNamesMappedFile != null ) {
            dumpQNames(locations.qNamesMappedFile, pipelineOptions, qNamesMultiMap);
        }
        log("Discovered " + qNamesMultiMap.size() + " mapped template names.");

        return new Tuple2<>(intervals, qNamesMultiMap);
    }

    private static void writeMetadata( final ReadMetadata readMetadata,
                                       final String filename,
                                       final PipelineOptions pipelineOptions ) {
        try ( final Writer writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(filename, pipelineOptions))) ) {
            writer.write("#reads:\t"+readMetadata.getNReads()+"\n");
            writer.write("#partitions:\t"+readMetadata.getNPartitions()+"\n");
            writer.write("max reads/partition:\t"+readMetadata.getMaxReadsInPartition()+"\n");
            writer.write("coverage:\t"+readMetadata.getCoverage()+"\n");
            for ( final Map.Entry<String, ReadMetadata.ReadGroupFragmentStatistics> entry :
                    readMetadata.getAllGroupStatistics().entrySet() ) {
                ReadMetadata.ReadGroupFragmentStatistics stats = entry.getValue();
                String name = entry.getKey();
                if ( name == null ) name = "NoGroup";
                writer.write("group "+name+":\t"+stats.getMedianFragmentSize()+
                        "-"+stats.getMedianNegativeDeviation()+"+"+stats.getMedianPositiveDeviation()+"\n");
            }
        }
        catch ( IOException ioe ) {
            throw new GATKException("Can't write metadata file.", ioe);
        }
    }

    /**
     * Kmerize each read mapped into a breakpoint interval,
     * get the template names of all reads sharing these kmers (regardless of where or if they're mapped),
     * and add these template names to the set of names for each interval.
     * Intervals having too many reads are killed.
     * The return is a description (as intervalId and explanatory String) of the intervals that were killed.
     */
    private Map<Integer, String> addAssemblyQNames(
            final Params params,
            final JavaSparkContext ctx,
            final String kmersToIgnoreFile,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> allPrimaryLines,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final JavaRDD<GATKRead> goodPrimaryLines =
                allPrimaryLines.filter(read -> !read.isDuplicate() && !read.failsVendorQualityCheck());

        final Tuple2<Map<Integer, String>, HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerAndIntervalsSet(params, ctx, kmersToIgnoreFile, qNamesMultiMap, nIntervals,
                                        goodPrimaryLines, locations, pipelineOptions);
        qNamesMultiMap.addAll(
                getAssemblyQNames(
                        params,
                        ctx,
                        kmerIntervalsAndDispositions._2(),
                        goodPrimaryLines));

        if ( locations.qNamesAssemblyFile != null ) {
            dumpQNames(locations.qNamesAssemblyFile, pipelineOptions, qNamesMultiMap);
        }

        log("Discovered "+qNamesMultiMap.size()+" unique template names for assembly.");
        return kmerIntervalsAndDispositions._1();
    }

    /**
     * Kmerize reads having template names in a given set,
     * filter out low complexity kmers and kmers that appear too often in the genome to be helpful in localizing reads,
     * kill intervals that have too few surviving kmers.
     * The return is a Tuple2 in which
     * _1 describes the intervals that have been killed for having too few kmers (as a map from intervalId onto an explanatory string),
     * and _2 describes the good kmers that we want to use in local assemblies (as a multimap from kmer onto intervalId).
     */
    private Tuple2<Map<Integer, String>, HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> getKmerAndIntervalsSet(
            final Params params,
            final JavaSparkContext ctx,
            final String kmersToIgnoreFile,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> goodPrimaryLines,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final Set<SVKmer> kmerKillSet = SVUtils.readKmersFile(params.kSize, kmersToIgnoreFile, pipelineOptions);
        log("Ignoring " + kmerKillSet.size() + " genomically common kmers.");

        final Tuple2<Map<Integer, String>, List<KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerIntervals(params, ctx, qNamesMultiMap, nIntervals, kmerKillSet, goodPrimaryLines, locations, pipelineOptions);
        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap =
                new HopscotchUniqueMultiMap<>(kmerIntervalsAndDispositions._2());
        log("Discovered " + kmerMultiMap.size() + " kmers.");

        return new Tuple2<>(kmerIntervalsAndDispositions._1(), kmerMultiMap);
    }

    /**
     * Transform all the reads for a supplied set of template names in each inverval into FASTQ records
     * for each interval, and do something with the list of FASTQ records for each interval (like write it to a file).
     */
    @VisibleForTesting static Map<Integer, String> generateFastqs(final JavaSparkContext ctx,
                                       final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
                                       final JavaRDD<GATKRead> reads,
                                       final int nIntervals,
                                       final boolean includeMappingLocation,
                                       final org.apache.spark.api.java.function.Function<Tuple2<Integer, List<byte[]>>, Tuple2<Integer, String>> fastqHandler) {
        final Broadcast<HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> broadcastQNamesMultiMap =
                ctx.broadcast(qNamesMultiMap);
        final int nPartitions = reads.partitions().size();

        final Map<Integer, String> intervalDispositions =
            reads
                .mapPartitionsToPair(readItr ->
                        new ReadsForQNamesFinder(broadcastQNamesMultiMap.value(), nIntervals,
                                includeMappingLocation).call(readItr).iterator(), false)
                .combineByKey(x -> x,
                                FindBreakpointEvidenceSpark::combineLists,
                                FindBreakpointEvidenceSpark::combineLists,
                                new HashPartitioner(nPartitions), false, null)
                .map(fastqHandler)
                .collect()
                .stream()
                .collect(Collectors.toMap(Tuple2::_1, Tuple2::_2));

        broadcastQNamesMultiMap.destroy();

        return intervalDispositions;
    }

    /** Concatenate two lists. */
    private static List<byte[]> combineLists( final List<byte[]> list1, final List<byte[]> list2 ) {
        final List<byte[]> result = new ArrayList<>(list1.size() + list2.size());
        result.addAll(list1);
        result.addAll(list2);
        return result;
    }

    /** write a FASTQ file for an assembly */
    @VisibleForTesting static Tuple2<Integer, String> writeFastq( final Tuple2<Integer, List<byte[]>> intervalAndFastqs,
                                    final String outputDir, final int maxFastqSize ) {
        final List<byte[]> fastqsList = intervalAndFastqs._2;
        SVFastqUtils.sortFastqRecords(fastqsList);

        final int fastqSize = fastqsList.stream().mapToInt(fastqRec -> fastqRec.length).sum();
        final String disposition;
        if ( fastqSize > maxFastqSize ) disposition = "FASTQ not written -- too big (" + fastqSize + " bytes).";
        else {
            final String fileName = outputDir + "/assembly" + intervalAndFastqs._1() + ".fastq";
            SVFastqUtils.writeFastqFile(fileName, null, fastqsList);
            disposition = fileName;
        }
        return new Tuple2<>(intervalAndFastqs._1(), disposition);
    }

    /**
     * Grab template names for all reads that contain kmers associated with a given breakpoint.
     */
    @VisibleForTesting static List<QNameAndInterval> getAssemblyQNames(
            final Params params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap,
            final JavaRDD<GATKRead> reads ) {
        final Broadcast<HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> broadcastKmerMultiMap =
                ctx.broadcast(kmerMultiMap);

        final int kSize = params.kSize;
        final double minEntropy = params.minEntropy;
        final int maxQNamesPerKmer = params.maxQNamesPerKmer;
        final int kmerMapSize = params.assemblyKmerMapSize;
        final List<QNameAndInterval> qNames =
            reads
                .mapPartitionsToPair(readItr ->
                        new MapPartitioner<>(readItr,
                                new QNamesForKmersFinder(kSize, minEntropy, broadcastKmerMultiMap.value())).iterator(), false)
                .mapPartitions(pairItr ->
                        new KmerQNameToQNameIntervalMapper(broadcastKmerMultiMap.value(),
                                                            maxQNamesPerKmer,
                                                            kmerMapSize).call(pairItr).iterator())
                .collect();

        broadcastKmerMultiMap.destroy();

        return qNames;
    }

    /** find kmers for each interval */
    @VisibleForTesting static Tuple2<Map<Integer, String>, List<KmerAndInterval>> getKmerIntervals(
            final Params params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final Set<SVKmer> kmerKillSet,
            final JavaRDD<GATKRead> reads,
            final Locations locations,
            final PipelineOptions pipelineOptions ) {

        final Broadcast<Set<SVKmer>> broadcastKmerKillSet = ctx.broadcast(kmerKillSet);
        final Broadcast<HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> broadcastQNameAndIntervalsMultiMap =
                ctx.broadcast(qNamesMultiMap);

        // given a set of template names with interval IDs and a kill set of ubiquitous kmers,
        // produce a set of interesting kmers for each interval ID
        final int kmersPerPartitionGuess = params.cleanerKmersPerPartitionGuess;
        final int minKmers = params.cleanerMinKmerCount;
        final int maxKmers = params.cleanerMaxKmerCount;
        final int maxIntervals = params.cleanerMaxIntervals;
        final int kSize = params.kSize;
        final double minEntropy = params.minEntropy;
        final List<KmerAndInterval> kmerIntervals =
            reads
                .mapPartitionsToPair(readItr ->
                        new MapPartitioner<>(readItr,
                            new QNameKmerizer(broadcastQNameAndIntervalsMultiMap.value(),
                                            broadcastKmerKillSet.value(), kSize, minEntropy)).iterator(), false)
                .reduceByKey(Integer::sum)
                .mapPartitions(itr -> new KmerCleaner(itr, kmersPerPartitionGuess, minKmers, maxKmers, maxIntervals).iterator())
                .collect();

        broadcastQNameAndIntervalsMultiMap.destroy();
        broadcastKmerKillSet.destroy();

        final int[] intervalKmerCounts = new int[nIntervals];
        for ( final KmerAndInterval kmerAndInterval : kmerIntervals ) {
            intervalKmerCounts[kmerAndInterval.getIntervalId()] += 1;
        }
        final Set<Integer> intervalsToKill = new HashSet<>();
        final Map<Integer, String> intervalDispositions = new HashMap<>();
        for ( int idx = 0; idx != nIntervals; ++idx ) {
            if ( intervalKmerCounts[idx] < params.minKmersPerInterval ) {
                intervalsToKill.add(idx);
                intervalDispositions.put(idx, "FASTQ not written -- too few kmers");
            }
        }

        final Iterator<QNameAndInterval> itr = qNamesMultiMap.iterator();
        while ( itr.hasNext() ) {
            if ( intervalsToKill.contains(itr.next().getIntervalId()) ) itr.remove();
        }

        final List<KmerAndInterval> filteredKmerIntervals = kmerIntervals.stream()
                .filter(kmerAndInterval -> !intervalsToKill.contains(kmerAndInterval.getIntervalId()))
                .collect(SVUtils.arrayListCollector(kmerIntervals.size()));

        // record the kmers with their interval IDs
        if ( locations.kmerFile != null ) {
            try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                    BucketUtils.createFile(locations.kmerFile, pipelineOptions)))) {
                for (final KmerAndInterval kmerAndInterval : filteredKmerIntervals) {
                    writer.write(kmerAndInterval.toString(kSize) + " " + kmerAndInterval.getIntervalId() + "\n");
                }
            } catch (final IOException ioe) {
                throw new GATKException("Can't write kmer intervals file " + locations.kmerFile, ioe);
            }
        }

        return new Tuple2<>(intervalDispositions, filteredKmerIntervals);
    }

    private List<SVInterval> removeIntervalsNearGapsAndLog( final List<SVInterval> intervals,
                                                            final int minDistanceToGap,
                                                            final ReadMetadata readMetadata,
                                                            final String exclusionIntervalsFile,
                                                            final PipelineOptions pipelineOptions ) {
        final List<SVInterval> result = removeIntervalsNearGaps(intervals, minDistanceToGap,
                readMetadata.getContigNameMap(), exclusionIntervalsFile, pipelineOptions);
        final int nKilledIntervals = intervals.size() - result.size();
        log("Killed " + nKilledIntervals + " intervals that were near reference gaps.");
        return result;
    }

    /** remove intervals that are near gaps */
    @VisibleForTesting static List<SVInterval> removeIntervalsNearGaps( final List<SVInterval> intervals,
                                                                        final int minDistanceToGap,
                                                                        final Map<String, Integer> contigNameMap,
                                                                        final String exclusionIntervalsFile,
                                                                        final PipelineOptions pipelineOptions ) {
        if ( exclusionIntervalsFile == null ) return intervals;
        final SortedSet<SVInterval> gaps =
                new TreeSet<>(SVUtils.readIntervalsFile(exclusionIntervalsFile, pipelineOptions, contigNameMap));
        return intervals.stream()
                .filter(interval -> {
                    final SortedSet<SVInterval> headSet = gaps.headSet(interval);
                    if ( !headSet.isEmpty() ) {
                        final SVInterval gapInterval = headSet.last();
                        if ( gapInterval.gapLen(interval) < minDistanceToGap ) return false;
                    }
                    final SortedSet<SVInterval> tailSet = gaps.tailSet(interval);
                    if ( !tailSet.isEmpty() ) {
                        final SVInterval gapInterval = tailSet.first();
                        if ( interval.gapLen(gapInterval) < minDistanceToGap ) return false;
                    }
                    return true;
                })
                .collect(Collectors.toCollection(() -> new ArrayList<>(intervals.size())));
    }

    private List<SVInterval> removeHighCoverageIntervalsAndLog( final Params params,
                                                                final JavaSparkContext ctx,
                                                                final Broadcast<ReadMetadata> broadcastMetadata,
                                                                final List<SVInterval> intervals,
                                                                final JavaRDD<GATKRead> mappedReads ) {
        final List<SVInterval> result = removeHighCoverageIntervals(params, ctx, broadcastMetadata, intervals, mappedReads);
        final int nKilledIntervals = intervals.size() - result.size();
        log("Killed " + nKilledIntervals + " intervals that had >" + params.maxIntervalCoverage + "x coverage.");
        return result;
    }

    /** figure out the coverage for each interval and filter out those with ridiculously high coverage */
    @VisibleForTesting static List<SVInterval> removeHighCoverageIntervals(
            final Params params,
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> reads) {
        final Broadcast<List<SVInterval>> broadcastIntervals = ctx.broadcast(intervals);
        final Map<Integer, Integer> intervalCoverage =
                reads
                    .mapPartitionsToPair(readItr ->
                            new IntervalCoverageFinder(broadcastMetadata.value(),broadcastIntervals.value(),readItr).iterator())
                    .reduceByKey(Integer::sum)
                    .collectAsMap();
        final int maxCoverage = params.maxIntervalCoverage;
        final List<SVInterval> result =
                IntStream.range(0, intervals.size())
                .filter(idx -> intervalCoverage.getOrDefault(idx, 0) <= maxCoverage*intervals.get(idx).getLength())
                .mapToObj(intervals::get)
                .collect(Collectors.toList());
        broadcastIntervals.destroy();
        return result;
    }

    /** find template names for reads mapping to each interval */
    @VisibleForTesting static HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> getQNames(
            final Params params,
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> reads ) {
        final Broadcast<List<SVInterval>> broadcastIntervals = ctx.broadcast(intervals);
        final List<QNameAndInterval> qNameAndIntervalList =
                reads
                        .mapPartitions(readItr ->
                                new MapPartitioner<>(readItr,
                                        new QNameFinder(broadcastMetadata.value(),
                                                broadcastIntervals.value())).iterator(), false)
                        .collect();
        broadcastIntervals.destroy();

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap =
                new HopscotchUniqueMultiMap<>(params.assemblyToMappedSizeRatioGuess*qNameAndIntervalList.size());
        qNamesMultiMap.addAll(qNameAndIntervalList);
        return qNamesMultiMap;
    }

    /** write template names and interval IDs to a file. */
    private static void dumpQNames( final String qNameFile,
                                    final PipelineOptions pipelineOptions,
                                    final Collection<QNameAndInterval> qNames ) {
        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(qNameFile, pipelineOptions)))) {
            for (final QNameAndInterval qnameAndInterval : qNames) {
                writer.write(qnameAndInterval.toString() + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write qname intervals file " + qNameFile, ioe);
        }
    }

    /**
     * Identify funky reads that support a hypothesis of a breakpoint in the vicinity, group the reads,
     * and declare a breakpoint interval where there is sufficient density of evidence.
     */
    @VisibleForTesting static List<SVInterval> getIntervals(
            final Params params,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> reads,
            final Locations locations ) {
        // find all breakpoint evidence, then filter for pile-ups
        final int maxFragmentSize = broadcastMetadata.value().getMaxMedianFragmentSize();
        final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
        final int nContigs = contigs.size();
        final int minMapQ = params.minEvidenceMapQ;
        final int minMatchLen = params.minEvidenceMatchLength;
        final int minEvidenceCount = params.minEvidenceCount;
        final JavaRDD<BreakpointEvidence> evidenceRDD =
                reads
                    .filter(read ->
                            read.getMappingQuality() >= minMapQ &&
                            read.getCigarElements()
                                    .stream()
                                    .filter(ele -> ele.getOperator()==CigarOperator.MATCH_OR_MISMATCH ||
                                                    ele.getOperator()==CigarOperator.EQ ||
                                                    ele.getOperator()==CigarOperator.X)
                                    .mapToInt(CigarElement::getLength)
                                    .sum() >= minMatchLen)
                    .mapPartitions(readItr ->
                            new MapPartitioner<>(readItr, new ReadClassifier(broadcastMetadata.value())).iterator(), true)
                    .mapPartitions(evidenceItr ->
                            new MapPartitioner<>(evidenceItr, new BreakpointClusterer(minEvidenceCount, 2*maxFragmentSize)).iterator(), true)
                    .mapPartitions(evidenceItr ->
                            new MapPartitioner<>(evidenceItr,
                                    new WindowSorter(3*maxFragmentSize), new BreakpointEvidence(nContigs)).iterator(), true);

        // record the evidence
        if ( locations.evidenceDir != null ) {
            evidenceRDD.cache();
            evidenceRDD.saveAsTextFile(locations.evidenceDir);
        }

        // find discrete intervals that contain the breakpoint evidence
        final Iterator<SVInterval> intervalItr =
                evidenceRDD
                        .mapPartitions(evidenceItr ->
                                new MapPartitioner<>(evidenceItr,
                                        new EvidenceToIntervalMapper(maxFragmentSize),
                                        new BreakpointEvidence(nContigs)).iterator(), true)
                        .collect()
                        .iterator();

        // coalesce overlapping intervals (can happen at partition boundaries)
        final List<SVInterval> intervals = new ArrayList<>();
        if ( intervalItr.hasNext() ) {
            SVInterval prev = intervalItr.next();
            while ( intervalItr.hasNext() ) {
                final SVInterval next = intervalItr.next();
                if ( prev.isDisjointFrom(next) ) {
                    intervals.add(prev);
                    prev = next;
                } else {
                    prev = prev.join(next);
                }
            }
            intervals.add(prev);
        }

        if ( locations.evidenceDir != null ) evidenceRDD.unpersist();

        return intervals;
    }

    private void log( final String message ) {
        logger.info(message);
    }

    @VisibleForTesting static class Locations {
        public final String metadataFile;
        public final String evidenceDir;
        public final String intervalFile;
        public final String qNamesMappedFile;
        public final String kmerFile;
        public final String qNamesAssemblyFile;
        public final String exclusionIntervalsFile;

        public Locations( final String metadataFile, final String evidenceDir, final String intervalFile,
                          final String qNamesMappedFile, final String kmerFile, final String qNamesAssemblyFile,
                          final String exclusionIntervalsFile) {
            this.metadataFile = metadataFile;
            this.evidenceDir = evidenceDir;
            this.intervalFile = intervalFile;
            this.qNamesMappedFile = qNamesMappedFile;
            this.kmerFile = kmerFile;
            this.qNamesAssemblyFile = qNamesAssemblyFile;
            this.exclusionIntervalsFile = exclusionIntervalsFile;
        }
    }

    @VisibleForTesting static class Params {
        public final int kSize;
        public final double minEntropy;
        public final int minEvidenceMapQ;
        public final int minEvidenceMatchLength;
        public final int maxIntervalCoverage;
        public final int minEvidenceCount;
        public final int minKmersPerInterval;
        public final int cleanerMaxIntervals;
        public final int cleanerMinKmerCount;
        public final int cleanerMaxKmerCount;
        public final int cleanerKmersPerPartitionGuess;
        public final int maxQNamesPerKmer;
        public final int assemblyKmerMapSize;
        public final int assemblyToMappedSizeRatioGuess;
        public final int maxFASTQSize;
        public final int exclusionIntervalPadding;

        public Params() {
            kSize = SVConstants.KMER_SIZE;          // kmer size
            minEntropy = SVConstants.MIN_ENTROPY;   // minimum kmer entropy
            minEvidenceMapQ = 20;                   // minimum map quality for evidential reads
            minEvidenceMatchLength = 45;            // minimum match length
            maxIntervalCoverage = 1000;             // maximum coverage on breakpoint interval
            minEvidenceCount = 15;                  // minimum number of evidentiary reads in called cluster
            minKmersPerInterval = 20;               // minimum number of good kmers in a valid interval
            cleanerMaxIntervals = 3;                // KmerCleaner maximum number of intervals a localizing kmer can appear in
            cleanerMinKmerCount = 3;                // KmerCleaner min kmer count
            cleanerMaxKmerCount = 125;              // KmerCleaner max kmer count
            cleanerKmersPerPartitionGuess = 600000; // KmerCleaner guess for number of unique error-free kmers per partition
            maxQNamesPerKmer = 500;                 // maximum template names for an assembly kmer
            assemblyKmerMapSize = 250000;           // guess for unique, error-free, scrubbed kmers per assembly partition
            assemblyToMappedSizeRatioGuess = 7;     // guess for ratio of total reads in assembly to evidentiary reads in interval
            maxFASTQSize = 10000000;                // maximum FASTQ size
            exclusionIntervalPadding = 0;           // exclusion interval extra padding
        }

        public Params( final int kSize, final double minEntropy, final int minEvidenceMapQ,
                       final int minEvidenceMatchLength, final int maxIntervalCoverage, final int minEvidenceCount,
                       final int minKmersPerInterval, final int cleanerMaxIntervals, final int cleanerMinKmerCount,
                       final int cleanerMaxKmerCount, final int cleanerKmersPerPartitionGuess,
                       final int maxQNamesPerKmer, final int asemblyKmerMapSize, final int assemblyToMappedSizeRatioGuess,
                       final int maxFASTQSize, final int exclusionIntervalPadding ) {
            this.kSize = kSize;
            this.minEntropy = minEntropy;
            this.minEvidenceMapQ = minEvidenceMapQ;
            this.minEvidenceMatchLength = minEvidenceMatchLength;
            this.maxIntervalCoverage = maxIntervalCoverage;
            this.minEvidenceCount = minEvidenceCount;
            this.minKmersPerInterval = minKmersPerInterval;
            this.cleanerMaxIntervals = cleanerMaxIntervals;
            this.cleanerMinKmerCount = cleanerMinKmerCount;
            this.cleanerMaxKmerCount = cleanerMaxKmerCount;
            this.cleanerKmersPerPartitionGuess = cleanerKmersPerPartitionGuess;
            this.maxQNamesPerKmer = maxQNamesPerKmer;
            this.assemblyKmerMapSize = asemblyKmerMapSize;
            this.assemblyToMappedSizeRatioGuess = assemblyToMappedSizeRatioGuess;
            this.maxFASTQSize = maxFASTQSize;
            this.exclusionIntervalPadding = exclusionIntervalPadding;
        }
    }

    /**
     * A class that acts as a filter for breakpoint evidence.
     * It passes only that evidence that is part of a putative cluster.
     */
    private static final class BreakpointClusterer
            implements Function<BreakpointEvidence, Iterator<BreakpointEvidence>> {
        private final int minEvidenceCount;
        private final int staleEventDistance;
        private final SortedMap<BreakpointEvidence, Boolean> locMap = new TreeMap<>();
        private final List<Map.Entry<BreakpointEvidence, Boolean>> reportableEntries;
        private final Iterator<BreakpointEvidence> noEvidence = Collections.emptyIterator();
        private int currentContig = -1;

        BreakpointClusterer( final int minEvidenceCount, final int staleEventDistance ) {
            this.minEvidenceCount = minEvidenceCount;
            this.staleEventDistance = staleEventDistance;
            this.reportableEntries = new ArrayList<>(2*minEvidenceCount);
        }

        @Override
        public Iterator<BreakpointEvidence> apply( final BreakpointEvidence evidence ) {
            if ( evidence.getContigIndex() != currentContig ) {
                currentContig = evidence.getContigIndex();
                locMap.clear();
            }

            locMap.put(evidence, true);

            final int locusStart = evidence.getEventStartPosition();
            final int locusEnd = evidence.getContigEnd();
            final int staleEnd = locusStart - staleEventDistance;
            int evidenceCount = 0;
            reportableEntries.clear();
            final Iterator<Map.Entry<BreakpointEvidence, Boolean>> itr = locMap.entrySet().iterator();
            while ( itr.hasNext() ) {
                final Map.Entry<BreakpointEvidence, Boolean> entry = itr.next();
                final BreakpointEvidence evidence2 = entry.getKey();
                final int contigEnd = evidence2.getContigEnd();
                if ( contigEnd <= staleEnd ) itr.remove();
                else if ( evidence2.getEventStartPosition() >= locusEnd ) break;
                else if ( contigEnd > locusStart ) {
                    evidenceCount += 1;
                    if ( entry.getValue() ) reportableEntries.add(entry);
                }
            }

            if ( evidenceCount >= minEvidenceCount ) {
                return reportableEntries.stream()
                        .map(entry -> { entry.setValue(false); return entry.getKey(); })
                        .iterator();
            }
            return noEvidence;
        }
    }

    /**
     * Class to fully sort a stream of nearly sorted BreakpointEvidences.
     */
    private static final class WindowSorter
            implements Function<BreakpointEvidence, Iterator<BreakpointEvidence>> {
        private final SortedSet<BreakpointEvidence> recordSet = new TreeSet<>();
        private final List<BreakpointEvidence> reportableEvidence = new ArrayList<>();
        private final int windowSize;
        private int currentContig = -1;

        WindowSorter( final int windowSize ) {
            this.windowSize = windowSize;
        }

        @Override
        public Iterator<BreakpointEvidence> apply( final BreakpointEvidence evidence ) {
            reportableEvidence.clear();
            if ( evidence.getContigIndex() != currentContig ) {
                reportableEvidence.addAll(recordSet);
                recordSet.clear();
                currentContig = evidence.getContigIndex();
            } else {
                final int reportableEnd = evidence.getEventStartPosition() - windowSize;
                final Iterator<BreakpointEvidence> itr = recordSet.iterator();
                while ( itr.hasNext() ) {
                    final BreakpointEvidence evidence2 = itr.next();
                    if ( evidence2.getEventStartPosition() >= reportableEnd ) break;
                    reportableEvidence.add(evidence2);
                    itr.remove();
                }
            }
            recordSet.add(evidence);
            return reportableEvidence.iterator();
        }
    }

    /**
     * A class to examine a stream of BreakpointEvidence, and group it into Intervals.
     */
    private static final class EvidenceToIntervalMapper implements Function<BreakpointEvidence, Iterator<SVInterval>> {
        private final Iterator<SVInterval> noInterval = Collections.emptyIterator();
        private final int gapSize;
        private int contig = -1;
        private int start;
        private int end;

        EvidenceToIntervalMapper(final int gapSize ) {
            this.gapSize = gapSize;
        }

        @Override
        public Iterator<SVInterval> apply( final BreakpointEvidence evidence ) {
            Iterator<SVInterval> result = noInterval;
            if ( evidence.getContigIndex() != contig ) {
                if ( contig != -1 ) {
                    result = new SingletonIterator<>(new SVInterval(contig, start, end));
                }
                contig = evidence.getContigIndex();
                start = evidence.getEventStartPosition();
                end = evidence.getContigEnd();
            } else if ( evidence.getEventStartPosition() >= end+gapSize ) {
                result = new SingletonIterator<>(new SVInterval(contig, start, end));
                start = evidence.getEventStartPosition();
                end = evidence.getContigEnd();
            } else {
                end = Math.max(end, evidence.getContigEnd());
            }
            return result;
        }
    }

    /**
     * Class to find the coverage of the intervals.
     */
    private static final class IntervalCoverageFinder implements Iterable<Tuple2<Integer, Integer>> {
        private final int[] basesInInterval;

        IntervalCoverageFinder( final ReadMetadata metadata,
                                final List<SVInterval> intervals,
                                final Iterator<GATKRead> readItr ) {
            basesInInterval = new int[intervals.size()];
            int intervalsIndex = 0;
            while ( readItr.hasNext() ) {
                final GATKRead read = readItr.next();
                final int readContigId = metadata.getContigID(read.getContig());
                final int readStart = read.getUnclippedStart();
                final int intervalsSize = intervals.size();
                while (intervalsIndex < intervalsSize) {
                    final SVInterval interval = intervals.get(intervalsIndex);
                    if (interval.getContig() > readContigId) break;
                    if (interval.getContig() == readContigId && interval.getEnd() > read.getStart()) break;
                    intervalsIndex += 1;
                }
                if (intervalsIndex >= intervalsSize) break;
                final SVInterval indexedInterval = intervals.get(intervalsIndex);
                final SVInterval readInterval = new SVInterval(readContigId, readStart, read.getUnclippedEnd());
                basesInInterval[intervalsIndex] += indexedInterval.overlapLen(readInterval);
            }
        }

        @Override
        public Iterator<Tuple2<Integer, Integer>> iterator() {
            return IntStream
                    .range(0, basesInInterval.length)
                    .filter(idx -> basesInInterval[idx] > 0)
                    .mapToObj(idx -> new Tuple2<>(idx, basesInInterval[idx]))
                    .iterator();
        }
    }

    /**
     * A template name and an intervalId.
     */
    @DefaultSerializer(QNameAndInterval.Serializer.class)
    @VisibleForTesting static final class QNameAndInterval implements Map.Entry<String, Integer> {
        private final byte[] qName;
        private final int hashVal;
        private final int intervalId;

        QNameAndInterval( final String qName, final int intervalId ) {
            this.qName = qName.getBytes();
            this.hashVal = qName.hashCode() ^ (47*intervalId);
            this.intervalId = intervalId;
        }

        private QNameAndInterval( final Kryo kryo, final Input input ) {
            final int nameLen = input.readInt();
            qName = input.readBytes(nameLen);
            hashVal = input.readInt();
            intervalId = input.readInt();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(qName.length);
            output.writeBytes(qName);
            output.writeInt(hashVal);
            output.writeInt(intervalId);
        }

        @Override
        public String getKey() { return getQName(); }
        @Override
        public Integer getValue() { return intervalId; }
        @Override
        public Integer setValue( final Integer value ) {
            throw new UnsupportedOperationException("Can't set QNameAndInterval.intervalId");
        }

        public String getQName() { return new String(qName); }
        public int getIntervalId() { return intervalId; }

        @Override
        public int hashCode() { return hashVal; }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof QNameAndInterval && equals((QNameAndInterval) obj);
        }

        public boolean equals( final QNameAndInterval that ) {
            return this.intervalId == that.intervalId && Arrays.equals(this.qName, that.qName);
        }

        public String toString() { return new String(qName)+" "+intervalId; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<QNameAndInterval> {
            @Override
            public void write( final Kryo kryo, final Output output, final QNameAndInterval qNameAndInterval ) {
                qNameAndInterval.serialize(kryo, output);
            }

            @Override
            public QNameAndInterval read( final Kryo kryo, final Input input, final Class<QNameAndInterval> klass ) {
                return new QNameAndInterval(kryo, input);
            }
        }
    }

    /**
     * Class to find the template names associated with reads in specified intervals.
     */
    private static final class QNameFinder implements Function<GATKRead, Iterator<QNameAndInterval>> {
        private final ReadMetadata metadata;
        private final List<SVInterval> intervals;
        private final Iterator<QNameAndInterval> noName = Collections.emptyIterator();
        private int intervalsIndex = 0;

        QNameFinder( final ReadMetadata metadata,
                     final List<SVInterval> intervals ) {
            this.metadata = metadata;
            this.intervals = intervals;
        }

        @Override
        public Iterator<QNameAndInterval> apply( final GATKRead read ) {
            final int readContigId = metadata.getContigID(read.getContig());
            final int readStart = read.getUnclippedStart();
            final int intervalsSize = intervals.size();
            while ( intervalsIndex < intervalsSize ) {
                final SVInterval interval = intervals.get(intervalsIndex);
                if ( interval.getContig() > readContigId ) break;
                if ( interval.getContig() == readContigId && interval.getEnd() > read.getStart() ) break;
                intervalsIndex += 1;
            }
            if ( intervalsIndex >= intervalsSize ) return noName;
            final SVInterval indexedInterval = intervals.get(intervalsIndex);
            final SVInterval readInterval = new SVInterval(readContigId, readStart, read.getUnclippedEnd());
            if ( indexedInterval.isDisjointFrom(readInterval) ) return noName;
            return new SingletonIterator<>(new QNameAndInterval(read.getName(), intervalsIndex));
        }
    }

    /**
     * A <Kmer,IntervalId> pair.
     */
    @DefaultSerializer(KmerAndInterval.Serializer.class)
    @VisibleForTesting final static class KmerAndInterval extends SVKmer implements Map.Entry<SVKmer, Integer> {
        private final int intervalId;

        KmerAndInterval( final SVKmer kmer, final int intervalId ) {
            super(kmer);
            this.intervalId = intervalId;
        }

        private KmerAndInterval( final Kryo kryo, final Input input ) {
            super(kryo, input);
            intervalId = input.readInt();
        }

        @Override
        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(intervalId);
        }

        @Override
        public SVKmer getKey() { return new SVKmer(this); }
        @Override
        public Integer getValue() { return intervalId; }
        @Override
        public Integer setValue( final Integer value ) {
            throw new UnsupportedOperationException("Can't set KmerAndInterval.intervalId");
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof KmerAndInterval && equals((KmerAndInterval)obj);
        }

        public boolean equals( final KmerAndInterval that ) {
            return super.equals(that) && this.intervalId == that.intervalId;
        }

        public int getIntervalId() { return intervalId; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<KmerAndInterval> {
            @Override
            public void write( final Kryo kryo, final Output output, final KmerAndInterval kmerAndInterval) {
                kmerAndInterval.serialize(kryo, output);
            }

            @Override
            public KmerAndInterval read(final Kryo kryo, final Input input,
                                        final Class<KmerAndInterval> klass ) {
                return new KmerAndInterval(kryo, input);
            }
        }
    }

    /**
     * Class that acts as a mapper from a stream of reads to a stream of KmerAndIntervals.
     * The template names of reads to kmerize, along with a set of kmers to ignore are passed in (by broadcast).
     */
    private static final class QNameKmerizer implements Function<GATKRead, Iterator<Tuple2<KmerAndInterval, Integer>>> {
        private final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNameAndIntervalMultiMap;
        private final Set<SVKmer> kmersToIgnore;
        private final int kSize;
        private final double minEntropy;
        private final ArrayList<Tuple2<KmerAndInterval, Integer>> tupleList = new ArrayList<>();

        QNameKmerizer( final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNameAndIntervalMultiMap,
                       final Set<SVKmer> kmersToIgnore, final int kSize, final double minEntropy ) {
            this.qNameAndIntervalMultiMap = qNameAndIntervalMultiMap;
            this.kmersToIgnore = kmersToIgnore;
            this.kSize = kSize;
            this.minEntropy = minEntropy;
        }

        @Override
        public Iterator<Tuple2<KmerAndInterval, Integer>> apply( final GATKRead read ) {
            final String qName = read.getName();
            final Iterator<QNameAndInterval> names = qNameAndIntervalMultiMap.findEach(qName);
            tupleList.clear();
            while ( names.hasNext() ) {
                final int intervalId = names.next().getIntervalId();
                SVKmerizerWithLowComplexityFilter.stream(read.getBases(), kSize, minEntropy)
                        .map(kmer -> kmer.canonical(kSize))
                        .filter(kmer -> !kmersToIgnore.contains(kmer))
                        .map(kmer -> new KmerAndInterval(kmer, intervalId))
                        .forEach(kmerCountAndInterval -> tupleList.add(new Tuple2<>(kmerCountAndInterval, 1)));
            }
            return tupleList.iterator();
        }
    }

    /**
     * Eliminates dups, and removes over-represented kmers.
     */
    private static final class KmerCleaner implements Iterable<KmerAndInterval> {

        private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap;

        KmerCleaner( final Iterator<Tuple2<KmerAndInterval, Integer>> kmerCountItr,
                     final int kmersPerPartitionGuess,
                     final int minKmerCount,
                     final int maxKmerCount,
                     final int maxIntervalsPerKmer ) {
            kmerMultiMap = new HopscotchUniqueMultiMap<>(kmersPerPartitionGuess);

            // remove kmers with extreme counts that won't help in building a local assembly
            while (kmerCountItr.hasNext()) {
                final Tuple2<KmerAndInterval, Integer> kmerCount = kmerCountItr.next();
                final int count = kmerCount._2;
                if (count >= minKmerCount && count <= maxKmerCount) kmerMultiMap.add(kmerCount._1);
            }

            final HopscotchSet<SVKmer> uniqueKmers = new HopscotchSet<>(kmerMultiMap.size());
            kmerMultiMap.stream().map(KmerAndInterval::getKey).forEach(uniqueKmers::add);
            uniqueKmers.stream()
                    .filter(kmer -> SVUtils.iteratorSize(kmerMultiMap.findEach(kmer)) > maxIntervalsPerKmer)
                    .forEach(kmerMultiMap::removeEach);
         }

        @Override
        public Iterator<KmerAndInterval> iterator() { return kmerMultiMap.iterator(); }
    }

    /**
     * Class that acts as a mapper from a stream of reads to a stream of <kmer,qname> pairs for a set of interesting kmers.
     * A multimap of interesting kmers is given to the constructor (by broadcast).
     */
    private static final class QNamesForKmersFinder implements Function<GATKRead, Iterator<Tuple2<SVKmer, String>>> {
        private final int kSize;
        private final double minEntropy;
        private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap;

        QNamesForKmersFinder( final int kSize, final double minEntropy,
                              final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap ) {
            this.kSize = kSize;
            this.minEntropy = minEntropy;
            this.kmerMultiMap = kmerMultiMap;
        }

        @Override
        public Iterator<Tuple2<SVKmer, String>> apply(final GATKRead read ) {
            List<Tuple2<SVKmer, String>> results = new ArrayList<>();
            SVKmerizerWithLowComplexityFilter.stream(read.getBases(), kSize, minEntropy)
                    .map( kmer -> kmer.canonical(kSize) )
                    .forEach( kmer -> {
                        final Iterator<KmerAndInterval> itr = kmerMultiMap.findEach(kmer);
                        if ( itr.hasNext() ) results.add(new Tuple2<>(kmer, read.getName()));
                    });
            return results.iterator();
        }
    }

    /**
     * Class that maps a stream of <kmer,qname> pairs into a stream of QNameAndIntervals.
     * A multimap of kmers onto intervalIds is given to the constructor.
     * Kmers that have too many (defined by constructor param) associated qnames are discarded.
     */
    private static final class KmerQNameToQNameIntervalMapper {
        private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap;
        private final int maxQNamesPerKmer;
        private final int kmerMapSize;

        KmerQNameToQNameIntervalMapper( final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap,
                                        final int maxQNamesPerKmer,
                                        final int kmerMapSize ) {
            this.kmerMultiMap = kmerMultiMap;
            this.maxQNamesPerKmer = maxQNamesPerKmer;
            this.kmerMapSize = kmerMapSize;
        }

        public Iterable<QNameAndInterval> call( final Iterator<Tuple2<SVKmer, String>> pairItr ) {
            HopscotchMap<SVKmer, List<String>, Map.Entry<SVKmer, List<String>>> kmerQNamesMap =
                    new HopscotchMap<>(kmerMapSize);
            while ( pairItr.hasNext() ) {
                final Tuple2<SVKmer, String> pair = pairItr.next();
                final SVKmer kmer = pair._1();
                Map.Entry<SVKmer, List<String>> entry = kmerQNamesMap.find(kmer);
                if ( entry == null ) {
                    // new entries are created with an empty list of qnames as their value,
                    // but if the list becomes too long we destroy it (by setting the value to null).
                    entry = new AbstractMap.SimpleEntry<>(kmer, new ArrayList<>());
                    kmerQNamesMap.add(entry);
                }
                List<String> qNames = entry.getValue();
                // if we're still growing the list
                if ( qNames != null ) {
                    // if the list becomes too long, discard it
                    if ( qNames.size() >= maxQNamesPerKmer ) entry.setValue(null);
                    else qNames.add(pair._2());
                }
            }

            final int qNameCount =
                    kmerQNamesMap.stream().mapToInt(entry -> entry.getValue()==null ? 0 : entry.getValue().size()).sum();
            HopscotchSet<QNameAndInterval> qNameAndIntervals = new HopscotchSet<>(qNameCount);
            for ( Map.Entry<SVKmer, List<String>> entry : kmerQNamesMap ) {
                final List<String> qNames = entry.getValue();
                // if the list hasn't been discarded for having grown too big
                if ( qNames != null ) {
                    Iterator<KmerAndInterval> intervalItr = kmerMultiMap.findEach(entry.getKey());
                    while ( intervalItr.hasNext() ) {
                        final int intervalId = intervalItr.next().getIntervalId();
                        for ( final String qName : qNames ) {
                            qNameAndIntervals.add(new QNameAndInterval(qName, intervalId));
                        }
                    }
                }
            }
            return qNameAndIntervals;
        }
    }

    /**
     * Find <intervalId,fastqBytes> pairs for interesting template names.
     */
    private static final class ReadsForQNamesFinder {
        private final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap;
        private final int nIntervals;
        private final int nReadsPerInterval;
        private final boolean includeMappingLocation;

        ReadsForQNamesFinder( final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
                              final int nIntervals, final boolean includeMappingLocation ) {
            this.qNamesMultiMap = qNamesMultiMap;
            this.nIntervals = nIntervals;
            this.nReadsPerInterval = 2*qNamesMultiMap.size()/nIntervals;
            this.includeMappingLocation = includeMappingLocation;
        }

        public Iterable<Tuple2<Integer, List<byte[]>>> call( final Iterator<GATKRead> readsItr ) {
            @SuppressWarnings({"unchecked", "rawtypes"})
            final List<byte[]>[] intervalReads = new List[nIntervals];
            int nPopulatedIntervals = 0;
            while ( readsItr.hasNext() ) {
                final GATKRead read = readsItr.next();
                final String readName = read.getName();
                final Iterator<QNameAndInterval> namesItr = qNamesMultiMap.findEach(readName);
                byte[] fastqBytes = null;
                while ( namesItr.hasNext() ) {
                    if ( fastqBytes == null ) fastqBytes = SVFastqUtils.readToFastqRecord(read, includeMappingLocation);
                    final int intervalId = namesItr.next().getIntervalId();
                    if ( intervalReads[intervalId] == null ) {
                        intervalReads[intervalId] = new ArrayList<>(nReadsPerInterval);
                        nPopulatedIntervals += 1;
                    }
                    intervalReads[intervalId].add(fastqBytes);
                }
            }
            final List<Tuple2<Integer, List<byte[]>>> fastQRecords = new ArrayList<>(nPopulatedIntervals);
            if ( nPopulatedIntervals > 0 ) {
                for ( int idx = 0; idx != nIntervals; ++idx ) {
                    final List<byte[]> readList = intervalReads[idx];
                    if ( readList != null ) fastQRecords.add(new Tuple2<>(idx, readList));
                }
            }
            return fastQRecords;
        }
    }
}
