package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.FlatMapGluer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexSingleton;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;
import static org.broadinstitute.hellbender.tools.spark.sv.evidence.BreakpointEvidence.ReadEvidence;

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
@BetaFeature
public final class FindBreakpointEvidenceSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;


    @ArgumentCollection
    private final FindBreakpointEvidenceSparkArgumentCollection params =
            new FindBreakpointEvidenceSparkArgumentCollection();


    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputAssemblyAlignments;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        gatherEvidenceAndWriteContigSamFile(ctx, params, getHeaderForReads(),
                getUnfilteredReads(), outputAssemblyAlignments, params.assembliesSortOrder, LogManager.getLogger(FindBreakpointEvidenceSpark.class));

    }

    /**
     * Gathers evidence reads and outputs them in a directory where reads are written out as interleaved FASTQ's.
     * Also produces the SAM records of contigs locally-assembled from such reads for downstream variant discovery.
     *
     * @return the in-memory representation of assembled contigs alignments, whose length equals the number of local assemblies (regardless of success of failure status)
     */
    public static List<AlignedAssemblyOrExcuse> gatherEvidenceAndWriteContigSamFile(
            final JavaSparkContext ctx,
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final String outputAssembliesFile,
            final SAMFileHeader.SortOrder outputSortOrder,
            final Logger toolLogger) {

        Utils.validate(header.getSortOrder() == SAMFileHeader.SortOrder.coordinate,
                "The reads must be coordinate sorted.");

        final SVReadFilter filter = new SVReadFilter(params);

        // develop evidence, intervals, and, finally, a set of template names for each interval
        final Tuple2<List<SVInterval>, HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> intervalsAndQNameMap =
                getMappedQNamesSet(params, ctx, header, unfilteredReads, filter, toolLogger);
        final List<SVInterval> intervals = intervalsAndQNameMap._1;
        if ( intervals.isEmpty() ) return new ArrayList<>();

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap = intervalsAndQNameMap._2;

        // supplement the template names with other reads that share kmers
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        if ( params.intervalOnlyAssembly ) {
            alignedAssemblyOrExcuseList = new ArrayList<>();
        } else {
            alignedAssemblyOrExcuseList = addAssemblyQNames(params, ctx, qNamesMultiMap, intervals.size(),
                    unfilteredReads, filter, toolLogger);
        }

        // write a FASTQ file for each interval
        final FermiLiteAssemblyHandler fermiLiteAssemblyHandler =
                new FermiLiteAssemblyHandler(params.alignerIndexImageFile, params.maxFASTQSize,
                                                params.fastqDir, params.gfaDir);
        alignedAssemblyOrExcuseList.addAll(
                handleAssemblies(ctx, qNamesMultiMap, unfilteredReads, filter, intervals.size(),
                        params.includeMappingLocation, fermiLiteAssemblyHandler));

        alignedAssemblyOrExcuseList.sort(Comparator.comparingInt(AlignedAssemblyOrExcuse::getAssemblyId));

        // record the intervals
        if ( params.intervalFile != null ) {
            AlignedAssemblyOrExcuse.writeIntervalFile(params.intervalFile, header, intervals, alignedAssemblyOrExcuseList);
        }

        // write the output file
        final SAMFileHeader cleanHeader = new SAMFileHeader(header.getSequenceDictionary());
        cleanHeader.setSortOrder(outputSortOrder);

        AlignedAssemblyOrExcuse.writeSAMFile(outputAssembliesFile, cleanHeader, alignedAssemblyOrExcuseList, outputSortOrder == SAMFileHeader.SortOrder.queryname);
        log("Wrote SAM file of aligned contigs.", toolLogger);

        return alignedAssemblyOrExcuseList;
    }

    /**
     * Find the breakpoint evidence,
     * cluster the evidence into intervals,
     * find the template names mapped into each interval,
     * kmerize each of these templates,
     * clean up by removing some intervals that are bogus as evidenced by ubiquitous kmers,
     * and return a set of template names and the intervals to which they belong.
     */
    private static Tuple2<List<SVInterval>, HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> getMappedQNamesSet(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Set<Integer> crossContigsToIgnoreSet;
        if ( params.crossContigsToIgnoreFile == null ) crossContigsToIgnoreSet = Collections.emptySet();
        else crossContigsToIgnoreSet = readCrossContigsToIgnoreFile(params.crossContigsToIgnoreFile,
                                                                    header.getSequenceDictionary());
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigsToIgnoreSet, header, params.maxTrackedFragmentLength, unfilteredReads, filter);
        if ( params.metadataFile != null ) {
            ReadMetadata.writeMetadata(readMetadata, params.metadataFile);
        }
        log("Metadata retrieved.", logger);

        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);
        List<SVInterval> intervals = getIntervals(params, broadcastMetadata, header, unfilteredReads, filter);

        final int nIntervals = intervals.size();
        log("Discovered " + nIntervals + " intervals.", logger);

        if ( nIntervals == 0 ) return new Tuple2<>(intervals, null);

        if ( params.exclusionIntervalsFile != null ) {
            intervals = removeIntervalsNearGapsAndLog(intervals, params.exclusionIntervalPadding, readMetadata,
                    params.exclusionIntervalsFile, logger);
        }
        intervals = removeHighCoverageIntervalsAndLog(
                                    params, ctx, broadcastMetadata, intervals, unfilteredReads, filter, logger);

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap =
                getQNames(params, ctx, broadcastMetadata, intervals, unfilteredReads, filter);
        broadcastMetadata.destroy();

        if ( params.qNamesMappedFile != null ) {
            QNameAndInterval.writeQNames(params.qNamesMappedFile, qNamesMultiMap);
        }
        log("Discovered " + qNamesMultiMap.size() + " mapped template names.", logger);

        return new Tuple2<>(intervals, qNamesMultiMap);
    }

    /** Read a file of contig names that will be ignored when checking for inter-contig pairs. */
    private static Set<Integer> readCrossContigsToIgnoreFile( final String crossContigsToIgnoreFile,
                                                              final SAMSequenceDictionary dictionary ) {
        final Set<Integer> ignoreSet = new HashSet<>();
        try ( final BufferedReader rdr =
                      new BufferedReader(
                              new InputStreamReader(BucketUtils.openFile(crossContigsToIgnoreFile))) ) {
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
    private static List<AlignedAssemblyOrExcuse> addAssemblyQNames(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Tuple2<List<AlignedAssemblyOrExcuse>, HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerAndIntervalsSet(params, ctx, qNamesMultiMap, nIntervals,
                                        unfilteredReads, filter, logger);
        qNamesMultiMap.addAll(
                getAssemblyQNames(params, ctx, kmerIntervalsAndDispositions._2(), unfilteredReads, filter));

        if ( params.qNamesAssemblyFile != null ) {
            QNameAndInterval.writeQNames(params.qNamesAssemblyFile, qNamesMultiMap);
        }

        log("Discovered "+qNamesMultiMap.size()+" unique template names for assembly.", logger);
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
    private static Tuple2<List<AlignedAssemblyOrExcuse>, HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> getKmerAndIntervalsSet(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Set<SVKmer> kmerKillSet =
                SVUtils.readKmersFile(params.kSize,
                                        params.kmersToIgnoreFile,
                                        new SVKmerLong(params.kSize));
        log("Ignoring " + kmerKillSet.size() + " genomically common kmers.", logger);

        final Tuple2<List<AlignedAssemblyOrExcuse>, List<KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerIntervals(params, ctx, qNamesMultiMap, nIntervals, kmerKillSet, unfilteredReads, filter);
        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap =
                new HopscotchUniqueMultiMap<>(kmerIntervalsAndDispositions._2());
        log("Discovered " + kmerMultiMap.size() + " kmers.", logger);

        return new Tuple2<>(kmerIntervalsAndDispositions._1(), kmerMultiMap);
    }

    /**
     * Functional interface that consumes the raw materials of an assembly to be aligned (i.e., a Tuple2 of assemblyId
     * and list of sequences) and returns an aligned assembly or an excuse for not producing one.
     */
    @VisibleForTesting interface LocalAssemblyHandler
            extends Serializable, Function<Tuple2<Integer,List<SVFastqUtils.FastqRead>>, AlignedAssemblyOrExcuse> {
    }

    /**
     * Transform all the reads for a supplied set of template names in each interval into FASTQ records
     * for each interval, and do something with the list of FASTQ records for each interval (like write it to a file).
     */
    @VisibleForTesting static List<AlignedAssemblyOrExcuse> handleAssemblies(
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final int nIntervals,
            final boolean includeMappingLocation,
            final LocalAssemblyHandler localAssemblyHandler ) {
        final Broadcast<HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> broadcastQNamesMultiMap =
                ctx.broadcast(qNamesMultiMap);
        final List<AlignedAssemblyOrExcuse> intervalDispositions =
            unfilteredReads
                .mapPartitionsToPair(readItr ->
                        new ReadsForQNamesFinder(broadcastQNamesMultiMap.value(), nIntervals,
                                includeMappingLocation, readItr, filter).iterator(), false)
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
    private static List<SVFastqUtils.FastqRead> combineLists(final List<SVFastqUtils.FastqRead> list1,
                                                             final List<SVFastqUtils.FastqRead> list2 ) {
        final List<SVFastqUtils.FastqRead> result = new ArrayList<>(list1.size() + list2.size());
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
        public AlignedAssemblyOrExcuse apply( final Tuple2<Integer, List<SVFastqUtils.FastqRead>> intervalAndReads ) {
            final List<SVFastqUtils.FastqRead> readsList = intervalAndReads._2();

            final int fastqSize = readsList.stream().mapToInt(FastqRead -> FastqRead.getBases().length).sum();
            if (fastqSize > maxFastqSize) {
                return new AlignedAssemblyOrExcuse(intervalAndReads._1(),
                        "no assembly -- too big (" + fastqSize + " bytes).");
            }
            if ( fastqDir != null ) {
                final String fastqName = String.format("%s/%s.fastq",fastqDir, AlignedAssemblyOrExcuse.formatAssemblyID(intervalAndReads._1()));
                final ArrayList<SVFastqUtils.FastqRead> sortedReads = new ArrayList<>(intervalAndReads._2());
                sortedReads.sort(Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
                SVFastqUtils.writeFastqFile(fastqName, sortedReads.iterator());
            }
            final FermiLiteAssembly assembly = new FermiLiteAssembler().createAssembly(readsList);
            if ( gfaDir != null ) {
                final String gfaName =  String.format("%s/%s.gfa",gfaDir, AlignedAssemblyOrExcuse.formatAssemblyID(intervalAndReads._1()));
                try ( final OutputStream os = BucketUtils.createFile(gfaName) ) {
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
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerMultiMap,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter ) {
        final Broadcast<HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> broadcastKmerMultiMap =
                ctx.broadcast(kmerMultiMap);

        final int kSize = params.kSize;
        final int maxDUSTScore = params.maxDUSTScore;
        final int maxQNamesPerKmer = params.maxQNamesPerKmer;
        final int kmerMapSize = params.assemblyKmerMapSize;
        final List<QNameAndInterval> qNames =
            unfilteredReads
                .mapPartitionsToPair(readItr ->
                        new FlatMapGluer<>(
                                new QNamesForKmersFinder(kSize, maxDUSTScore, broadcastKmerMultiMap.value(), filter),
                                readItr), false)
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
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final Set<SVKmer> kmerKillSet,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter ) {

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
            unfilteredReads
                .mapPartitionsToPair(readItr ->
                        new FlatMapGluer<>(
                                new QNameKmerizer(
                                        broadcastQNameAndIntervalsMultiMap.value(),
                                        broadcastKmerKillSet.value(), kSize, maxDUSTScore, filter),
                                readItr), false)
                .reduceByKey(Integer::sum)
                .mapPartitions(itr ->
                        new KmerCleaner(itr, kmersPerPartitionGuess, minKmers, maxKmers, maxIntervals).iterator())
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
        if ( params.kmerFile != null ) {
            try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                    BucketUtils.createFile(params.kmerFile)))) {
                for (final KmerAndInterval kmerAndInterval : filteredKmerIntervals) {
                    writer.write(kmerAndInterval.toString(kSize) + " " + kmerAndInterval.getIntervalId() + "\n");
                }
            } catch (final IOException ioe) {
                throw new GATKException("Can't write kmer intervals file " + params.kmerFile, ioe);
            }
        }

        return new Tuple2<>(intervalDispositions, filteredKmerIntervals);
    }

    private static List<SVInterval> removeIntervalsNearGapsAndLog( final List<SVInterval> intervals,
                                                                   final int minDistanceToGap,
                                                                   final ReadMetadata readMetadata,
                                                                   final String exclusionIntervalsFile,
                                                                   final Logger logger ) {
        final List<SVInterval> result = removeIntervalsNearGaps(intervals, minDistanceToGap,
                readMetadata.getContigNameMap(), exclusionIntervalsFile);
        final int nKilledIntervals = intervals.size() - result.size();
        log("Killed " + nKilledIntervals + " intervals that were near reference gaps.", logger);
        return result;
    }

    /** remove intervals that are near gaps */
    @VisibleForTesting static List<SVInterval> removeIntervalsNearGaps( final List<SVInterval> intervals,
                                                                        final int minDistanceToGap,
                                                                        final Map<String, Integer> contigNameMap,
                                                                        final String exclusionIntervalsFile ) {
        if ( exclusionIntervalsFile == null ) return intervals;
        final SortedSet<SVInterval> gaps =
                new TreeSet<>(SVUtils.readIntervalsFile(exclusionIntervalsFile, contigNameMap));
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

    private static List<SVInterval> removeHighCoverageIntervalsAndLog(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger) {
        final List<SVInterval> result =
                removeHighCoverageIntervals(params, ctx, broadcastMetadata, intervals, unfilteredReads, filter);
        final int nKilledIntervals = intervals.size() - result.size();
        log("Killed " + nKilledIntervals + " intervals that had >" + params.maxIntervalCoverage + "x coverage.", logger);
        return result;
    }

    /** figure out the coverage for each interval and filter out those with ridiculously high coverage */
    @VisibleForTesting static List<SVInterval> removeHighCoverageIntervals(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter ) {
        final Broadcast<List<SVInterval>> broadcastIntervals = ctx.broadcast(intervals);
        final Map<Integer, Integer> intervalCoverage =
                unfilteredReads
                    .mapPartitionsToPair(readItr ->
                        new IntervalCoverageFinder(
                                broadcastMetadata.value(),broadcastIntervals.value(),readItr,filter).iterator())
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
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<SVInterval> intervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter ) {
        final Broadcast<List<SVInterval>> broadcastIntervals = ctx.broadcast(intervals);
        final List<QNameAndInterval> qNameAndIntervalList =
                unfilteredReads
                    .mapPartitions(readItr ->
                            new FlatMapGluer<>(
                                    new QNameFinder(broadcastMetadata.value(), broadcastIntervals.value(), filter),
                                    readItr), false)
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
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter ) {
        // find all breakpoint evidence, then filter for pile-ups
        final int nContigs = header.getSequenceDictionary().getSequences().size();
        final int minEvidenceWeight = params.minEvidenceWeight;
        final int minCoherentEvidenceWeight = params.minCoherentEvidenceWeight;
        final int allowedOverhang = params.allowedShortFragmentOverhang;

        // 1) identify well-mapped reads
        // 2) that look like they support a hypothesis of a breakpoint in the vicinity
        // 3a) filter out those that lack supporting evidence from a sufficient number of other reads, except
        // 3b) pass through everything within a fragment length of partition boundaries
        final JavaRDD<BreakpointEvidence> evidenceRDD =
                unfilteredReads
                    .mapPartitions(readItr -> {
                        final GATKRead sentinel = new SAMRecordToGATKReadAdapter(null);
                        return FlatMapGluer.applyMapFunc(
                                new ReadClassifier(broadcastMetadata.value(),sentinel,allowedOverhang,filter),
                                readItr,sentinel);
                        }, true)
                    .mapPartitionsWithIndex( (idx, evidenceItr) -> {
                            final ReadMetadata readMetadata = broadcastMetadata.value();
                            final PartitionCrossingChecker xChecker =
                                    new PartitionCrossingChecker(idx, readMetadata,
                                                                 readMetadata.getMaxMedianFragmentSize());
                            return new BreakpointDensityFilter(evidenceItr,readMetadata,
                                                               minEvidenceWeight,minCoherentEvidenceWeight,xChecker);
                        }, true);

        evidenceRDD.cache();

        // record the evidence
        if ( params.evidenceDir != null ) {
            evidenceRDD
                    .filter(BreakpointEvidence::isValidated)
                    .saveAsTextFile(params.evidenceDir);
        }

        // replace clumps of evidence with a single aggregate indicator of evidence (to save memory), except
        // don't clump anything within two fragment lengths of a partition boundary.
        // collect the whole mess in the driver.
        final int maxFragmentSize = broadcastMetadata.value().getMaxMedianFragmentSize();
        final List<BreakpointEvidence> collectedEvidence =
                evidenceRDD
                        .mapPartitionsWithIndex( (idx, readEvidenceItr) ->
                                new FlatMapGluer<>(
                                        new BreakpointEvidenceClusterer(maxFragmentSize,
                                                new PartitionCrossingChecker(idx,broadcastMetadata.value(),2*maxFragmentSize)),
                                        readEvidenceItr,
                                        new BreakpointEvidence(new SVInterval(nContigs,1,1),0,false)), true)
                        .collect();

        evidenceRDD.unpersist();

        // reapply the density filter (all data collected -- no more worry about partition boundaries).
        final Iterator<BreakpointEvidence> evidenceIterator =
                new BreakpointDensityFilter(collectedEvidence.iterator(),
                        broadcastMetadata.value(), minEvidenceWeight, minCoherentEvidenceWeight,
                        new PartitionCrossingChecker());
        final List<BreakpointEvidence> allEvidence = new ArrayList<>(collectedEvidence.size());
        while ( evidenceIterator.hasNext() ) {
            allEvidence.add(evidenceIterator.next());
        }

        // write additional validated read evidence from partition-boundary reads
        if ( params.evidenceDir != null ) {
            final String crossPartitionFile = params.evidenceDir+"/part-xxxxx";
            try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(crossPartitionFile))) ) {
                for ( final BreakpointEvidence ev : allEvidence ) {
                    // Only tell 'em about the ReadEvidence that was validated in the driver's pass over the stream.
                    // (The validated ReadEvidence instances well away from the partition boundaries that have already
                    // been reported will have been replaced by a generic BreakpointEvidence object that represents
                    // a whole cluster.  A few bits of validated ReadEvidence near the partition boundaries will
                    // get double reported, unfortunately.)
                    if ( ev instanceof ReadEvidence ) {
                        writer.write(ev.toString());
                        writer.write('\n');
                    }
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write cross-partition evidence to "+crossPartitionFile, ioe);
            }
        }

        // re-cluster the new evidence -- we can now glue across partition boundaries
        Iterator<BreakpointEvidence> evidenceIterator2 =
            new FlatMapGluer<>(new BreakpointEvidenceClusterer(maxFragmentSize,new PartitionCrossingChecker()),
                               allEvidence.iterator(),
                               new BreakpointEvidence(new SVInterval(nContigs,1,1),0,false));

        final List<SVInterval> intervals = new ArrayList<>(allEvidence.size());
        while ( evidenceIterator2.hasNext() ) {
            intervals.add(evidenceIterator2.next().getLocation());
        }

        return intervals;
    }

    private static void log(final String message, final Logger logger) {
        logger.info(message);
    }
}
