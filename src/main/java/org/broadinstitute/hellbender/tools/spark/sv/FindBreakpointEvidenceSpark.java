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
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
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
        programGroup = SparkProgramGroup.class)
public final class FindBreakpointEvidenceSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    //--------- parameters ----------
    @VisibleForTesting static final Params defaultParams =
            new Params(
                    SVConstants.KMER_SIZE,   // kmer size
                    SVConstants.MIN_ENTROPY, // minimum kmer entropy
                    20,                      // minimum map quality for evidential reads
                    45,                      // minimum match length
                    1000,                    // maximum coverage on breakpoint interval
                    15,                      // minimum number of evidentiary reads in called cluster
                    20000000,                // guess for total kmers in BAM partition
                    2,                       // minimum kmer count within partition of high frequency kmer
                    250000,                  // guess for number of unique error-free kmers in one partition
                    200,                     // high frequency kmer count
                    3,                       // maximum number of intervals a localizing kmer can appear in
                    3,                       // KmerCleaner min kmer count
                    125,                     // KmerCleaner max kmer count
                    600000,                  // KmerCleaner guess for number of unique error-free kmers per partition
                    7,                       // guess for ratio of total reads in assembly to evidentiary reads in interval
                    100000000,               // maximum FASTQ size
                    0                        // exclusion interval extra padding
            );

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

    @Argument(doc = "Guess for the total number of kmers in one partition of the input file.",
            fullName = "totalKmersPerPartitionGuess")
    private int totalKmersPerPartitionGuess = defaultParams.totalKmersPerPartitionGuess;

    @Argument(doc = "Minimum count of kmer within reads partition to be considered in finding high frequency kmers.",
            fullName = "minKmerCountWithinPartition")
    private int minKmerCountWithinPartition = defaultParams.minKmerCountWithinPartition;

    @Argument(doc = "Unique error-free kmers per partition", fullName = "uniqueErrorFreeKmersPerPartition")
    private int uniqueErrorFreeKmersPerPartitionGuess = defaultParams.uniqueErrorFreeKmersPerPartitionGuess;

    @Argument(doc = "Count for kmer to be considered high frequency.", fullName = "minHighFrequencyKmerCount")
    private int minHighFrequencyKmerCount = defaultParams.minHighFrequencyKmerCount;

    @Argument(doc = "KmerCleaner maximum number of intervals for a localizing kmer.", fullName = "cleanerMaxIntervals")
    private int cleanerMaxIntervals = defaultParams.cleanerMaxIntervals;

    @Argument(doc = "KmerCleaner minimum kmer count.", fullName = "cleanerMinKmerCount")
    private int cleanerMinKmerCount = defaultParams.cleanerMinKmerCount;

    @Argument(doc = "KmerCleaner maximum kmer count.", fullName = "cleanerMaxKmerCount")
    private int cleanerMaxKmerCount = defaultParams.cleanerMaxKmerCount;

    @Argument(doc = "KmerCleaner unique error-free kmers per partition", fullName = "cleanerKmersPerPartitionGuess")
    private int cleanerKmersPerPartitionGuess = defaultParams.cleanerKmersPerPartitionGuess;

    @Argument(doc = "Guess at the ratio of reads in the final assembly to the number reads mapped to the interval.",
            fullName = "assemblyToMappedSizeRatioGuess")
    private int assemblyToMappedSizeRatioGuess = defaultParams.assemblyToMappedSizeRatioGuess;

    @Argument(doc = "Maximum FASTQ file size.", fullName = "maxFASTQSize")
    private int maxFASTQSize = defaultParams.maxFASTQSize;

    @Argument(doc = "Exclusion interval padding.", fullName = "exclusionIntervalPadding")
    private int exclusionIntervalPadding = defaultParams.exclusionIntervalPadding;

    @Argument(doc = "Include read mapping location in FASTQ files.", fullName = "includeMappingLocation")
    private boolean includeMappingLocation = true;

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

    @Argument(doc = "file for high frequency kmers output", fullName = "highFrequencyKmers", optional = true)
    private String highFrequencyKmersFile;

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
                                highFrequencyKmersFile, kmerFile, qNamesAssemblyFile, exclusionIntervalsFile);
        final Params params =
                new Params(kSize, minEntropy, minEvidenceMapQ, minEvidenceMatchLength, maxIntervalCoverage,
                            minEvidenceCount, totalKmersPerPartitionGuess, minKmerCountWithinPartition,
                            uniqueErrorFreeKmersPerPartitionGuess, minHighFrequencyKmerCount,
                            cleanerMaxIntervals, cleanerMinKmerCount, cleanerMaxKmerCount, cleanerKmersPerPartitionGuess,
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
        addAssemblyQNames(params, ctx, kmersToIgnoreFile, qNamesMultiMap, allPrimaryLines, locations, pipelineOptions);

        // write a FASTQ file for each interval
        final String outDir = outputDir;
        final int maxFastqSize = maxFASTQSize;
        final boolean includeMapLoc = includeMappingLocation;
        final List<Tuple2<Integer, String>> intervalDispositions =
                generateFastqs(ctx, qNamesMultiMap, allPrimaryLines, intervals.size(), includeMapLoc,
                                intervalAndFastqBytes -> writeFastq(intervalAndFastqBytes, outDir, maxFastqSize));

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
                                           final List<Tuple2<Integer, String>> intervalDispositions ) {
        final int nIntervals = intervals.size();
        final String[] dispositions = new String[nIntervals];
        for ( final Tuple2<Integer, String> tuple : intervalDispositions ) {
            final int intervalId = tuple._1();
            if ( intervalId < 0 || intervalId >= nIntervals ) {
                throw new GATKException("Unexpected intervalId among dispositions: "+intervalId);
            }
            dispositions[tuple._1()] = tuple._2();
        }

        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(intervalFile, pipelineOptions)))) {
            final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
            for ( int intervalId = 0; intervalId != nIntervals; ++intervalId ) {
                final SVInterval interval = intervals.get(intervalId);
                final String seqName = contigs.get(interval.getContig()).getSequenceName();
                String disposition = dispositions[intervalId];
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
     */
    private void addAssemblyQNames(
            final Params params,
            final JavaSparkContext ctx,
            final String kmersToIgnoreFile,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final JavaRDD<GATKRead> allPrimaryLines,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final JavaRDD<GATKRead> goodPrimaryLines =
                allPrimaryLines.filter(read -> !read.isDuplicate() && !read.failsVendorQualityCheck());

        qNamesMultiMap.addAll(
                getAssemblyQNames(
                        params,
                        ctx,
                        getKmerAndIntervalsSet(params, ctx, kmersToIgnoreFile, qNamesMultiMap,
                                                goodPrimaryLines, locations, pipelineOptions),
                        goodPrimaryLines));

        if ( locations.qNamesAssemblyFile != null ) {
            dumpQNames(locations.qNamesAssemblyFile, pipelineOptions, qNamesMultiMap);
        }

        log("Discovered "+qNamesMultiMap.size()+" unique template names for assembly.");
    }

    /**
     * Kmerize reads having template names in a given set,
     * filter out kmers that appear too often in this read set or in the genome to be helpful in localizing reads,
     * and return the set of kmers that appear in each interval.
     */
    private HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> getKmerAndIntervalsSet(
            final Params params,
            final JavaSparkContext ctx,
            final String kmersToIgnoreFile,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final JavaRDD<GATKRead> goodPrimaryLines,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final Set<SVKmer> kmerKillSet = SVUtils.readKmersFile(params.kSize, kmersToIgnoreFile, pipelineOptions);
        log("Ignoring " + kmerKillSet.size() + " genomically common kmers.");
        final List<SVKmer> kmerKillList = getHighCountKmers(params, goodPrimaryLines, locations, pipelineOptions);
        log("Ignoring " + kmerKillList.size() + " common kmers in the reads.");
        kmerKillSet.addAll(kmerKillList);
        log("Ignoring a total of " + kmerKillSet.size() + " unique common kmers.");

        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap =
                new HopscotchUniqueMultiMap<>(
                    getKmerIntervals(params, ctx, qNamesMultiMap, kmerKillSet, goodPrimaryLines, locations, pipelineOptions));
        log("Discovered " + kmerMultiMap.size() + " kmers.");

        return kmerMultiMap;
    }

    /**
     * Transform all the reads for a supplied set of template names in each inverval into FASTQ records
     * for each interval, and do something with the list of FASTQ records for each interval (like write it to a file).
     */
    @VisibleForTesting static List<Tuple2<Integer, String>> generateFastqs(final JavaSparkContext ctx,
                                       final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
                                       final JavaRDD<GATKRead> reads,
                                       final int nIntervals,
                                       final boolean includeMappingLocation,
                                       final org.apache.spark.api.java.function.Function<Tuple2<Integer, List<byte[]>>, Tuple2<Integer, String>> fastqHandler) {
        final Broadcast<HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> broadcastQNamesMultiMap =
                ctx.broadcast(qNamesMultiMap);
        final int nPartitions = reads.partitions().size();

        final List<Tuple2<Integer, String>> intervalDispositions =
            reads
                .mapPartitionsToPair(readItr ->
                        new ReadsForQNamesFinder(broadcastQNamesMultiMap.value(), nIntervals,
                                includeMappingLocation).call(readItr), false)
                .combineByKey(x -> x,
                                FindBreakpointEvidenceSpark::combineLists,
                                FindBreakpointEvidenceSpark::combineLists,
                                new HashPartitioner(nPartitions), false, null)
                .map(fastqHandler)
                .collect();

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
        final List<QNameAndInterval> qNames =
            reads
                .mapPartitions(readItr ->
                        new MapPartitioner<>(readItr,
                                new QNamesForKmersFinder(kSize, minEntropy, broadcastKmerMultiMap.value())), false)
                .collect();

        broadcastKmerMultiMap.destroy();

        return qNames;
    }

    /** find kmers for each interval */
    @VisibleForTesting static List<KmerAndInterval> getKmerIntervals(
            final Params params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
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
                                            broadcastKmerKillSet.value(), kSize, minEntropy)), false)
                .reduceByKey(Integer::sum)
                .mapPartitions(itr -> new KmerCleaner(itr, kmersPerPartitionGuess, minKmers, maxKmers, maxIntervals))
                .collect();

        broadcastQNameAndIntervalsMultiMap.destroy();
        broadcastKmerKillSet.destroy();

        // record the kmers with their interval IDs
        if ( locations.kmerFile != null ) {
            try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                    BucketUtils.createFile(locations.kmerFile, pipelineOptions)))) {
                for (final KmerAndInterval kmerAndInterval : kmerIntervals) {
                    writer.write(kmerAndInterval.toString(kSize) + " " + kmerAndInterval.getIntervalId() + "\n");
                }
            } catch (final IOException ioe) {
                throw new GATKException("Can't write kmer intervals file " + locations.kmerFile, ioe);
            }
        }

        return kmerIntervals;
    }

    /** Find all the kmers in this read set that occur with high frequency. */
    @VisibleForTesting static List<SVKmer> getHighCountKmers(
            final Params params,
            final JavaRDD<GATKRead> reads,
            final Locations locations,
            final PipelineOptions pipelineOptions ) {
        final int nPartitions = reads.partitions().size();
        final int kSize = params.kSize;
        final double minEntropy = params.minEntropy;
        final int totalKmersPerPartitionGuess = params.totalKmersPerPartitionGuess;
        final int minKmerCountWithinPartition = params.minKmerCountWithinPartition;
        final int uniqueErrorFreeKmersPerPartitionGuess = params.uniqueErrorFreeKmersPerPartitionGuess;
        final int minHighFrequencyKmerCount = params.minHighFrequencyKmerCount;
        final List<SVKmer> kmers =
                reads
                    .mapPartitions(readItr ->
                            new KmerCounter(readItr, kSize, minEntropy, totalKmersPerPartitionGuess, minKmerCountWithinPartition), false)
                    .mapToPair(kmerAndCount -> new Tuple2<>(kmerAndCount, null))
                    .partitionBy(new HashPartitioner(nPartitions))
                    .map(tuple -> tuple._1)
                    .mapPartitions(kmerItr ->
                            new KmerReducer(kmerItr, uniqueErrorFreeKmersPerPartitionGuess, minHighFrequencyKmerCount))
                    .map(SVKmer::new)
                    .collect();

        if ( locations.highFrequencyKmersFile != null ) {
            SVUtils.writeKmersFile(params.kSize, locations.highFrequencyKmersFile, pipelineOptions, kmers);
        }

        return kmers;
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
                            new IntervalCoverageFinder(broadcastMetadata.value(),broadcastIntervals.value(),readItr))
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
                                                broadcastIntervals.value())), false)
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
                            new MapPartitioner<>(readItr, new ReadClassifier(broadcastMetadata.value())), true)
                    .mapPartitions(evidenceItr ->
                            new MapPartitioner<>(evidenceItr, new BreakpointClusterer(minEvidenceCount, 2*maxFragmentSize)), true)
                    .mapPartitions(evidenceItr ->
                            new MapPartitioner<>(evidenceItr,
                                    new WindowSorter(3*maxFragmentSize), new BreakpointEvidence(nContigs)), true);

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
                                        new BreakpointEvidence(nContigs)), true)
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
        public final String highFrequencyKmersFile;
        public final String kmerFile;
        public final String qNamesAssemblyFile;
        public final String exclusionIntervalsFile;

        public Locations( final String metadataFile, final String evidenceDir, final String intervalFile,
                          final String qNamesMappedFile, final String highFrequencyKmersFile,
                          final String kmerFile, final String qNamesAssemblyFile,
                          final String exclusionIntervalsFile) {
            this.metadataFile = metadataFile;
            this.evidenceDir = evidenceDir;
            this.intervalFile = intervalFile;
            this.qNamesMappedFile = qNamesMappedFile;
            this.highFrequencyKmersFile = highFrequencyKmersFile;
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
        public final int totalKmersPerPartitionGuess;
        public final int minKmerCountWithinPartition;
        public final int uniqueErrorFreeKmersPerPartitionGuess;
        public final int minHighFrequencyKmerCount;
        public final int cleanerMaxIntervals;
        public final int cleanerMinKmerCount;
        public final int cleanerMaxKmerCount;
        public final int cleanerKmersPerPartitionGuess;
        public final int assemblyToMappedSizeRatioGuess;
        public final int maxFASTQSize;
        public final int exclusionIntervalPadding;

        public Params( final int kSize, final double minEntropy, final int minEvidenceMapQ,
                       final int minEvidenceMatchLength, final int maxIntervalCoverage, final int minEvidenceCount,
                       final int totalKmersPerPartitionGuess, final int minKmerCountWithinPartition,
                       final int uniqueErrorFreeKmersPerPartitionGuess, final int minHighFrequencyKmerCount,
                       final int cleanerMaxIntervals, final int cleanerMinKmerCount,
                       final int cleanerMaxKmerCount, final int cleanerKmersPerPartitionGuess,
                       final int assemblyToMappedSizeRatioGuess, final int maxFASTQSize,
                       final int exclusionIntervalPadding ) {
            this.kSize = kSize;
            this.minEntropy = minEntropy;
            this.minEvidenceMapQ = minEvidenceMapQ;
            this.minEvidenceMatchLength = minEvidenceMatchLength;
            this.maxIntervalCoverage = maxIntervalCoverage;
            this.minEvidenceCount = minEvidenceCount;
            this.totalKmersPerPartitionGuess = totalKmersPerPartitionGuess;
            this.minKmerCountWithinPartition = minKmerCountWithinPartition;
            this.uniqueErrorFreeKmersPerPartitionGuess = uniqueErrorFreeKmersPerPartitionGuess;
            this.minHighFrequencyKmerCount = minHighFrequencyKmerCount;
            this.cleanerMaxIntervals = cleanerMaxIntervals;
            this.cleanerMinKmerCount = cleanerMinKmerCount;
            this.cleanerMaxKmerCount = cleanerMaxKmerCount;
            this.cleanerKmersPerPartitionGuess = cleanerKmersPerPartitionGuess;
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

        public Iterator<SVInterval> apply(final BreakpointEvidence evidence ) {
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

        KmerAndInterval(final SVKmer kmer, final int intervalId ) {
            super(kmer);
            this.intervalId = intervalId;
        }

        private KmerAndInterval(final Kryo kryo, final Input input ) {
            super(kryo, input);
            intervalId = input.readInt();
        }

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
     * A <Kmer,count> pair.
     */
    @DefaultSerializer(KmerAndCount.Serializer.class)
    private final static class KmerAndCount extends SVKmer implements Map.Entry<SVKmer, Integer> {
        private int count;

        KmerAndCount( final SVKmer kmer ) {
            super(kmer);
            this.count = 1;
        }

        private KmerAndCount(final Kryo kryo, final Input input ) {
            super(kryo, input);
            count = input.readInt();
        }

        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(count);
        }

        @Override
        public SVKmer getKey() { return new SVKmer(this); }
        @Override
        public Integer getValue() { return count; }
        @Override
        public Integer setValue( final Integer value ) { final int oldCount = count; count = value; return oldCount; }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof KmerAndCount && equals((KmerAndCount)obj);
        }

        public boolean equals( final KmerAndCount that ) {
            return super.equals(that);
        }

        public int getCount() { return count; }
        public void incrementCount() { count += 1; }
        public void incrementCount( final int val ) { count += val; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<KmerAndCount> {
            @Override
            public void write( final Kryo kryo, final Output output, final KmerAndCount kmerAndCount) {
                kmerAndCount.serialize(kryo, output);
            }

            @Override
            public KmerAndCount read(final Kryo kryo, final Input input, final Class<KmerAndCount> klass ) {
                return new KmerAndCount(kryo, input);
            }
        }
    }

    /**
     * Kmerizes reads into a big Set, ignores those that occur only once (i.e., < MIN_KMER_COUNT),
     * and returns a Kmer and occurrence count for the rest.
     */
    private final static class KmerCounter implements Iterable<KmerAndCount> {
        private final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerMap;

        KmerCounter( final Iterator<GATKRead> readItr,
                     final int kSize,
                     final double minEntropy,
                     final int totalKmersPerPartitionGuess,
                     final int minKmerCountWithinPartition ) {
            final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerMap = new HopscotchMap<>(totalKmersPerPartitionGuess);
            while ( readItr.hasNext() ) {
                SVKmerizerWithLowComplexityFilter.stream(readItr.next().getBases(), kSize, minEntropy)
                        .map(kmer -> kmer.canonical(kSize))
                        .forEach(kmer -> {
                            final KmerAndCount kmerAndCount = kmerMap.find(kmer);
                            if ( kmerAndCount != null ) kmerAndCount.incrementCount();
                            else kmerMap.add(new KmerAndCount(kmer));
                        });
            }
            final Iterator<KmerAndCount> kmerItr = kmerMap.iterator();
            while ( kmerItr.hasNext() ) {
                if ( kmerItr.next().getCount() < minKmerCountWithinPartition ) kmerItr.remove();
            }
            this.kmerMap = kmerMap;
        }

        @Override
        public Iterator<KmerAndCount> iterator() { return kmerMap.iterator(); }
    }

    /**
     * Returns kmers that occur very frequently (>= MAX_KMER_COUNT).
     */
    private final static class KmerReducer implements Iterable<KmerAndCount> {
        private final HopscotchMap<SVKmer, Integer, KmerAndCount> kmerMap;

        KmerReducer( final Iterator<KmerAndCount> kmerItr,
                     final int uniqueErrorFreeKmersPerPartitionGuess, final int minHighFrequencyKmerCount ) {
            kmerMap = new HopscotchMap<>(uniqueErrorFreeKmersPerPartitionGuess);
            while ( kmerItr.hasNext() ) {
                final KmerAndCount kmerAndCount = kmerItr.next();
                final KmerAndCount tableKmerAndCount = kmerMap.find(kmerAndCount.getKey());
                if ( tableKmerAndCount != null ) tableKmerAndCount.incrementCount(kmerAndCount.getCount());
                else kmerMap.add(kmerAndCount);
            }
            final Iterator<KmerAndCount> kmerItr2 = kmerMap.iterator();
            while ( kmerItr2.hasNext() ) {
                if ( kmerItr2.next().getCount() < minHighFrequencyKmerCount ) kmerItr2.remove();
            }
        }

        @Override
        public Iterator<KmerAndCount> iterator() { return kmerMap.iterator(); }
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

        public Iterator<KmerAndInterval> iterator() { return kmerMultiMap.iterator(); }
    }

    /**
     * Class that acts as a mapper from a stream of reads to a stream of <intervalId,read> pairs.
     * It knows which breakpoint(s) a read belongs to (if any) by kmerizing the read, and looking up each SVKmer in
     * a multi-map of SVKmers onto intervalId.
     */
    private static final class QNamesForKmersFinder implements Function<GATKRead, Iterator<QNameAndInterval>> {
        private final int kSize;
        private final double minEntropy;
        private final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap;
        private final Set<Integer> intervalIdSet = new HashSet<>();
        private final List<QNameAndInterval> qNameAndIntervalList = new ArrayList<>();
        private final Iterator<QNameAndInterval> emptyIterator = Collections.emptyIterator();

        QNamesForKmersFinder( final int kSize, final double minEntropy,
                              final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap ) {
            this.kSize = kSize;
            this.minEntropy = minEntropy;
            this.kmerMultiMap = kmerMultiMap;
        }

        public Iterator<QNameAndInterval> apply(final GATKRead read) {
            intervalIdSet.clear();
            SVKmerizerWithLowComplexityFilter.stream(read.getBases(), kSize, minEntropy)
                    .map( kmer -> kmer.canonical(kSize) )
                    .forEach( kmer -> {
                        final Iterator<KmerAndInterval> itr = kmerMultiMap.findEach(kmer);
                        while ( itr.hasNext() ) {
                            intervalIdSet.add(itr.next().getValue());
                        }
                    });
            if (intervalIdSet.isEmpty()) return emptyIterator;
            qNameAndIntervalList.clear();
            final String qName = read.getName();
            intervalIdSet.forEach(intervalId ->
                    qNameAndIntervalList.add(new QNameAndInterval(qName, intervalId)));
            return qNameAndIntervalList.iterator();
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

        @SuppressWarnings("unchecked")
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
