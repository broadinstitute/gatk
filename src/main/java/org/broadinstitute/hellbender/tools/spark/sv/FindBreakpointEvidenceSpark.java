package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.SVFastqUtils.FastqRead;
import org.broadinstitute.hellbender.tools.spark.utils.*;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tool to discover reads that support a hypothesis of a genomic breakpoint.
 * Reads sharing kmers with reads aligned near putative breakpoints are pulled out
 *  for local assemblies of these breakpoint regions.
 * The local assemblies are done with FermiLite, and the assembled contigs are aligned to reference.
 * Final output is a SAM file of aligned contigs to be called for structural variants.
 */
@CommandLineProgramProperties(summary="Find reads that evidence breakpoints."+
        "  Pull reads for local assemblies in breakpoint regions using shared kmers."+
        "  Assemble breakpoint regions with FermiLite, and align assembled contigs to reference.",
        oneLineSummary="Prepare local assemblies of putative genomic breakpoints for structural variant discovery.",
        programGroup=StructuralVariationSparkProgramGroup.class)
public final class FindBreakpointEvidenceSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    //--------- parameters ----------

    // no-arg constructor for Params object establishes default values
    @VisibleForTesting static final Params defaultParams = new Params();

    @Argument(doc = "Kmer size.", fullName = "kSize")
    private int kSize = defaultParams.kSize;

    @Argument(doc = "maximum kmer DUST score", fullName = "kmerMaxDUSTScore")
    private int maxDUSTScore = SVConstants.MAX_DUST_SCORE;

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

    @Argument(doc = "Don't look for extra reads mapped outside the interval.", fullName = "intervalOnlyAssembly")
    private boolean intervalOnlyAssembly = false;

    // --------- locations ----------

    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputSAM;

    @Argument(doc = "bwa-mem index image file", fullName = "alignerIndexImage")
    private String alignerIndexImageFile;

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

    @Argument(doc = "output dir for assembled fastqs", fullName = "fastqDir", optional = true)
    private String fastqDir;

    @Argument(doc = "output dir for assemblies", fullName = "gfaDir", optional = true)
    private String gfaDir;

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

    /**
     * This is a path to a text file of contig names (one per line) that will be ignored when looking for inter-contig pairs.
     */
    @Argument(doc = "file containing alt contig names that will be ignored when looking for inter-contig pairs",
            fullName = "crossContigsToIgnore", optional = true)
    private String crossContigsToIgnoreFile;

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
                                kmerFile, qNamesAssemblyFile, exclusionIntervalsFile, alignerIndexImageFile);
        final Params params =
                new Params(kSize, maxDUSTScore, minEvidenceMapQ, minEvidenceMatchLength, maxIntervalCoverage,
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
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        if ( intervalOnlyAssembly ) {
            alignedAssemblyOrExcuseList = new ArrayList<>();
        } else {
            alignedAssemblyOrExcuseList = addAssemblyQNames(params, ctx, kmersToIgnoreFile, qNamesMultiMap, intervals.size(),
                    allPrimaryLines, locations, pipelineOptions);
        }

        // write a FASTQ file for each interval
        final FermiLiteAssemblyHandler fermiLiteAssemblyHandler =
                new FermiLiteAssemblyHandler(locations.alignerIndexImageFile, maxFASTQSize, fastqDir, gfaDir);
        alignedAssemblyOrExcuseList.addAll(
                handleAssemblies(ctx, qNamesMultiMap, allPrimaryLines, intervals.size(),
                                includeMappingLocation, fastqDir != null, fermiLiteAssemblyHandler));

        alignedAssemblyOrExcuseList.sort(Comparator.comparingInt(AlignedAssemblyOrExcuse::getAssemblyId));

        // record the intervals
        if ( locations.intervalFile != null ) {
            AlignedAssemblyOrExcuse.writeIntervalFile(locations.intervalFile, pipelineOptions, header, intervals, alignedAssemblyOrExcuseList);
        }

        // write the output file
        final SAMFileHeader cleanHeader = new SAMFileHeader(header.getSequenceDictionary());
        AlignedAssemblyOrExcuse.writeSAMFile(outputSAM, pipelineOptions, cleanHeader, alignedAssemblyOrExcuseList);
        log("Wrote SAM file of aligned contigs.");
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
        final Set<Integer> crossContigsToIgnoreSet;
        if ( crossContigsToIgnoreFile == null ) crossContigsToIgnoreSet = Collections.emptySet();
        else crossContigsToIgnoreSet = readCrossContigsToIgnoreFile(crossContigsToIgnoreFile,
                                                                    pipelineOptions,
                                                                    header.getSequenceDictionary());
        final JavaRDD<GATKRead> mappedReads =
                unfilteredReads.filter(read ->
                        !read.isDuplicate() && !read.failsVendorQualityCheck() && !read.isUnmapped());
        final ReadMetadata readMetadata = new ReadMetadata(crossContigsToIgnoreSet, header, mappedReads);
        if ( locations.metadataFile != null ) {
            ReadMetadata.writeMetadata(readMetadata, locations.metadataFile, pipelineOptions);
        }
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
            QNameAndInterval.writeQNames(locations.qNamesMappedFile, pipelineOptions, qNamesMultiMap);
        }
        log("Discovered " + qNamesMultiMap.size() + " mapped template names.");

        return new Tuple2<>(intervals, qNamesMultiMap);
    }

    /** Read a file of contig names that will be ignored when checking for inter-contig pairs. */
    private static Set<Integer> readCrossContigsToIgnoreFile( final String crossContigsToIgnoreFile,
                                                              final PipelineOptions pipelineOptions,
                                                              final SAMSequenceDictionary dictionary ) {
        final Set<Integer> ignoreSet = new HashSet<>();
        try ( final BufferedReader rdr =
                      new BufferedReader(
                              new InputStreamReader(BucketUtils.openFile(crossContigsToIgnoreFile,pipelineOptions))) ) {
            String line;
            while ( (line = rdr.readLine()) != null ) {
                final int tigId = dictionary.getSequenceIndex(line);
                if ( tigId == -1 ) {
                    throw new UserException("crossContigToIgnoreFile contains an unrecognized contig name: "+line);
                }
                ignoreSet.add(tigId);
            }
        }
        catch ( final IOException ioe ) {
            throw new UserException("Can't read crossContigToIgnore file "+crossContigsToIgnoreFile, ioe);
        }
        return ignoreSet;
    }

    /**
     * Kmerize each read mapped into a breakpoint interval,
     * get the template names of all reads sharing these kmers (regardless of where or if they're mapped),
     * and add these template names to the set of names for each interval.
     * Intervals having too many reads are killed.
     * The return is a description (as intervalId and explanatory String) of the intervals that were killed.
     */
    private List<AlignedAssemblyOrExcuse> addAssemblyQNames(
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

        final Tuple2<List<AlignedAssemblyOrExcuse>, HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerAndIntervalsSet(params, ctx, kmersToIgnoreFile, qNamesMultiMap, nIntervals,
                                        goodPrimaryLines, locations, pipelineOptions);
        qNamesMultiMap.addAll(
                getAssemblyQNames(
                        params,
                        ctx,
                        kmerIntervalsAndDispositions._2(),
                        goodPrimaryLines));

        if ( locations.qNamesAssemblyFile != null ) {
            QNameAndInterval.writeQNames(locations.qNamesAssemblyFile, pipelineOptions, qNamesMultiMap);
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
    private Tuple2<List<AlignedAssemblyOrExcuse>, HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> getKmerAndIntervalsSet(
            final Params params,
            final JavaSparkContext ctx,
            final String kmersToIgnoreFile,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> goodPrimaryLines,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final Set<SVKmer> kmerKillSet = SVUtils.readKmersFile(params.kSize, kmersToIgnoreFile, pipelineOptions, new SVKmerLong(params.kSize));
        log("Ignoring " + kmerKillSet.size() + " genomically common kmers.");

        final Tuple2<List<AlignedAssemblyOrExcuse>, List<KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerIntervals(params, ctx, qNamesMultiMap, nIntervals, kmerKillSet, goodPrimaryLines, locations, pipelineOptions);
        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap =
                new HopscotchUniqueMultiMap<>(kmerIntervalsAndDispositions._2());
        log("Discovered " + kmerMultiMap.size() + " kmers.");

        return new Tuple2<>(kmerIntervalsAndDispositions._1(), kmerMultiMap);
    }

    /**
     * Functional interface that consumes the raw materials of an assembly to be aligned (i.e., a Tuple2 of assemblyId
     * and list of sequences) and returns an aligned assembly or an excuse for not producing one.
     */
    @VisibleForTesting interface LocalAssemblyHandler
            extends Serializable, Function<Tuple2<Integer,List<FastqRead>>, AlignedAssemblyOrExcuse> {
    }

    /**
     * Transform all the reads for a supplied set of template names in each interval into FASTQ records
     * for each interval, and do something with the list of FASTQ records for each interval (like write it to a file).
     */
    @VisibleForTesting static List<AlignedAssemblyOrExcuse> handleAssemblies(
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final JavaRDD<GATKRead> reads,
            final int nIntervals,
            final boolean includeMappingLocation,
            final boolean dumpFASTQs,
            final LocalAssemblyHandler localAssemblyHandler ) {
        final Broadcast<HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> broadcastQNamesMultiMap =
                ctx.broadcast(qNamesMultiMap);
        final List<AlignedAssemblyOrExcuse> intervalDispositions =
            reads
                .mapPartitionsToPair(readItr ->
                        new ReadsForQNamesFinder(broadcastQNamesMultiMap.value(), nIntervals,
                                includeMappingLocation, dumpFASTQs).call(readItr).iterator(), false)
                .combineByKey(x -> x,
                                FindBreakpointEvidenceSpark::combineLists,
                                FindBreakpointEvidenceSpark::combineLists,
                                new HashPartitioner(nIntervals), false, null)
                .map(localAssemblyHandler::apply)
                .collect();

        broadcastQNamesMultiMap.destroy();
        BwaMemIndexSingleton.closeAllDistributedInstances(ctx);

        return intervalDispositions;
    }

    /** Concatenate two lists. */
    private static List<FastqRead> combineLists( final List<FastqRead> list1, final List<FastqRead> list2 ) {
        final List<FastqRead> result = new ArrayList<>(list1.size() + list2.size());
        result.addAll(list1);
        result.addAll(list2);
        return result;
    }

    /** This LocalAssemblyHandler aligns assembly contigs with BWA, along with some optional writing of intermediate results. */
    private static final class FermiLiteAssemblyHandler implements LocalAssemblyHandler {
        private static final long serialVersionUID = 1L;
        private final String alignerIndexFile;
        private final int maxFastqSize;
        private final String fastqDir;
        private final String gfaDir;

        FermiLiteAssemblyHandler( final String alignerIndexFile, final int maxFastqSize,
                                  final String fastqDir, final String gfaDir ) {
            this.alignerIndexFile = alignerIndexFile;
            this.maxFastqSize = maxFastqSize;
            this.fastqDir = fastqDir;
            this.gfaDir = gfaDir;
        }

        @Override
        public AlignedAssemblyOrExcuse apply( final Tuple2<Integer, List<FastqRead>> intervalAndReads ) {
            final List<FastqRead> readsList = intervalAndReads._2();

            final int fastqSize = readsList.stream().mapToInt(FastqRead -> FastqRead.getBases().length).sum();
            if (fastqSize > maxFastqSize) {
                return new AlignedAssemblyOrExcuse(intervalAndReads._1(),
                        "no assembly -- too big (" + fastqSize + " bytes).");
            }
            final PipelineOptions pipelineOptions = null;
            if ( fastqDir != null ) {
                final String fastqName = String.format("%s/assembly%06d.fastq",fastqDir,intervalAndReads._1());
                final ArrayList<FastqRead> sortedReads = new ArrayList<>(intervalAndReads._2());
                sortedReads.sort(Comparator.comparing(FastqRead::getName));
                SVFastqUtils.writeFastqFile(fastqName, pipelineOptions, sortedReads.iterator());
            }
            final FermiLiteAssembly assembly = new FermiLiteAssembler().createAssembly(readsList);
            if ( gfaDir != null ) {
                final String gfaName = String.format("%s/assembly%06d.gfa",gfaDir,intervalAndReads._1());
                try ( final OutputStream os = BucketUtils.createFile(gfaName, pipelineOptions) ) {
                    assembly.writeGFA(os);
                }
                catch ( final IOException ioe ) {
                    throw new GATKException("Can't write "+gfaName, ioe);
                }
            }
            final List<byte[]> tigSeqs =
                    assembly.getContigs().stream()
                            .map(FermiLiteAssembly.Contig::getSequence)
                            .collect(SVUtils.arrayListCollector(assembly.getNContigs()));
            try ( final BwaMemAligner aligner = new BwaMemAligner(BwaMemIndexSingleton.getInstance(alignerIndexFile)) ) {
                aligner.setIntraCtgOptions();
                final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(tigSeqs);
                return new AlignedAssemblyOrExcuse(intervalAndReads._1(), assembly, alignments);
            }
        }
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
        final int maxDUSTScore = params.maxDUSTScore;
        final int maxQNamesPerKmer = params.maxQNamesPerKmer;
        final int kmerMapSize = params.assemblyKmerMapSize;
        final List<QNameAndInterval> qNames =
            reads
                .mapPartitionsToPair(readItr ->
                        new MapPartitioner<>(readItr,
                                new QNamesForKmersFinder(kSize, maxDUSTScore, broadcastKmerMultiMap.value())).iterator(), false)
                .mapPartitions(pairItr ->
                        new KmerQNameToQNameIntervalMapper(broadcastKmerMultiMap.value(),
                                                            maxQNamesPerKmer,
                                                            kmerMapSize).call(pairItr).iterator())
                .collect();

        broadcastKmerMultiMap.destroy();

        return qNames;
    }

    /** find kmers for each interval */
    @VisibleForTesting static Tuple2<List<AlignedAssemblyOrExcuse>, List<KmerAndInterval>> getKmerIntervals(
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
        final int maxDUSTScore = params.maxDUSTScore;
        final List<KmerAndInterval> kmerIntervals =
            reads
                .mapPartitionsToPair(readItr ->
                        new MapPartitioner<>(readItr,
                            new QNameKmerizer(broadcastQNameAndIntervalsMultiMap.value(),
                                            broadcastKmerKillSet.value(), kSize, maxDUSTScore)).iterator(), false)
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
        final List<AlignedAssemblyOrExcuse> intervalDispositions = new ArrayList<>();
        for ( int idx = 0; idx != nIntervals; ++idx ) {
            if ( intervalKmerCounts[idx] < params.minKmersPerInterval ) {
                intervalsToKill.add(idx);
                intervalDispositions.add(new AlignedAssemblyOrExcuse(idx, "FASTQ not written -- too few kmers"));
            }
        }

        qNamesMultiMap.removeIf( qNameAndInterval -> intervalsToKill.contains(qNameAndInterval.getIntervalId()) );

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
        public final String alignerIndexImageFile;

        public Locations( final String metadataFile, final String evidenceDir, final String intervalFile,
                          final String qNamesMappedFile, final String kmerFile, final String qNamesAssemblyFile,
                          final String exclusionIntervalsFile, final String alignerIndexImageFile ) {
            this.metadataFile = metadataFile;
            this.evidenceDir = evidenceDir;
            this.intervalFile = intervalFile;
            this.qNamesMappedFile = qNamesMappedFile;
            this.kmerFile = kmerFile;
            this.qNamesAssemblyFile = qNamesAssemblyFile;
            this.exclusionIntervalsFile = exclusionIntervalsFile;
            this.alignerIndexImageFile = alignerIndexImageFile;
        }
    }

    @VisibleForTesting static class Params {
        public final int kSize;
        public final int maxDUSTScore;
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
            maxDUSTScore = SVConstants.MAX_DUST_SCORE;// maximum for DUST-like kmer complexity score
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
            maxFASTQSize = 3000000;                 // maximum assembly size (total input bases)
            exclusionIntervalPadding = 0;           // exclusion interval extra padding
        }

        public Params( final int kSize, final int maxDUSTScore, final int minEvidenceMapQ,
                       final int minEvidenceMatchLength, final int maxIntervalCoverage, final int minEvidenceCount,
                       final int minKmersPerInterval, final int cleanerMaxIntervals, final int cleanerMinKmerCount,
                       final int cleanerMaxKmerCount, final int cleanerKmersPerPartitionGuess,
                       final int maxQNamesPerKmer, final int asemblyKmerMapSize, final int assemblyToMappedSizeRatioGuess,
                       final int maxFASTQSize, final int exclusionIntervalPadding ) {
            this.kSize = kSize;
            this.maxDUSTScore = maxDUSTScore;
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
}
