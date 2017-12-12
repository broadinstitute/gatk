package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.FlatMapGluer;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.spark.SparkUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;
import static org.broadinstitute.hellbender.tools.spark.sv.evidence.BreakpointEvidence.ExternalEvidence;
import static org.broadinstitute.hellbender.tools.spark.sv.evidence.BreakpointEvidence.ReadEvidence;

/**
 * Tool to discover genomic breakpoints associated with structural variants.
 *
 * <p>This tool prepares local assemblies of putative genomic breakpoints for structural variant discovery.
 * Reads sharing kmers with reads aligned near putative breakpoints are pulled out for local assemblies of
 * these breakpoint regions.</p>
 * <p>The local assemblies are done with FermiLite, and the assembled contigs are aligned to reference.
 * Output is a file of aligned contigs to be used in calling structural variants.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>A file of paired-end, aligned and coordinate-sorted reads.</li>
 *     <li>A BWA index image for the reference.</li>
 *     <li>A list of ubiquitous kmers to ignore.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A file of aligned contigs.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk FindBreakpointEvidenceSpark \
 *     -I input_reads.bam \
 *     --aligner-index-image reference.img \
 *     --kmers-to-ignore ignored_kmers.txt \
 *     -O assemblies.sam
 * </pre>
 *
 * <h3>Notes</h3>
 * <p>Expected input is a paired-end, coordinate-sorted BAM with around 30x coverage.
 * Coverage much lower than that probably won't work well.</p>
 * <p>You can use BwaMemIndexImageCreator to create the index image file, and FindBadGenomicGenomicKmersSpark
 * to create the list of kmers to ignore.</p>
 */
@DocumentedFeature
@BetaFeature
@CommandLineProgramProperties(
        oneLineSummary = "Tool to discover genomic breakpoints associated with structural variants.",
        summary =
        "This tool prepares local assemblies of putative genomic breakpoints for structural variant discovery." +
        " Reads sharing kmers with reads aligned near putative breakpoints are pulled out for local assemblies of" +
        " these breakpoint regions." +
        " The local assemblies are done with FermiLite, and the assembled contigs are aligned to reference." +
        " Output is a file of aligned contigs to be used in calling structural variants.",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
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

        gatherEvidenceAndWriteContigSamFile(ctx, params, getHeaderForReads(), getUnfilteredReads(),
                outputAssemblyAlignments, logger);

    }

    /**
     * Gathers evidence reads and outputs them in a directory where reads are written out as interleaved FASTQ's.
     * Also produces the SAM records of contigs locally-assembled from such reads for downstream variant discovery.
     *
     * @return the in-memory representation of assembled contigs alignments, whose length equals the number of local assemblies (regardless of success of failure status)
     */
    public static AssembledEvidenceResults gatherEvidenceAndWriteContigSamFile(
            final JavaSparkContext ctx,
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final String outputAssemblyAlignments,
            final Logger logger) {

        Utils.validate(header.getSortOrder() == SAMFileHeader.SortOrder.coordinate,
                "The reads must be coordinate sorted.");

        final SVReadFilter filter = new SVReadFilter(params);

        final Set<Integer> crossContigsToIgnoreSet;
        if ( params.crossContigsToIgnoreFile == null ) crossContigsToIgnoreSet = Collections.emptySet();
        else crossContigsToIgnoreSet = readCrossContigsToIgnoreFile(params.crossContigsToIgnoreFile,
                header.getSequenceDictionary());
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigsToIgnoreSet, header, params.maxTrackedFragmentLength, unfilteredReads, filter, logger);
        if ( params.metadataFile != null ) {
            ReadMetadata.writeMetadata(readMetadata, params.metadataFile);
        }
        log("Metadata retrieved.", logger);

        // develop evidence, intervals, and, finally, a set of template names for each interval
        final EvidenceScanResults
                evidenceScanResults = getMappedQNamesSet(params, readMetadata, ctx, header, unfilteredReads, filter, logger);
        final List<SVInterval> intervals = evidenceScanResults.intervals;
        if ( intervals.isEmpty() ) return new AssembledEvidenceResults(
                evidenceScanResults.readMetadata,
                intervals,
                new ArrayList<>(),
                evidenceScanResults.evidenceTargetLinks);

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap = evidenceScanResults.qNamesForAssemblyMultiMap;

        // supplement the template names with other reads that share kmers
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        if ( params.intervalOnlyAssembly ) {
            alignedAssemblyOrExcuseList = new ArrayList<>();
        } else {
            alignedAssemblyOrExcuseList = addAssemblyQNames(params, readMetadata, ctx, qNamesMultiMap, intervals.size(),
                    unfilteredReads, filter, logger);
        }

        // write a FASTQ file for each interval
        final FermiLiteAssemblyHandler fermiLiteAssemblyHandler =
                new FermiLiteAssemblyHandler(params.alignerIndexImageFile, params.maxFASTQSize, params.fastqDir, params.writeGFAs);
        alignedAssemblyOrExcuseList.addAll(
                handleAssemblies(ctx, qNamesMultiMap, unfilteredReads, filter, intervals.size(),
                        params.includeMappingLocation, fermiLiteAssemblyHandler));

        alignedAssemblyOrExcuseList.sort(Comparator.comparingInt(AlignedAssemblyOrExcuse::getAssemblyId));

        // record the intervals
        if ( params.intervalFile != null ) {
            AlignedAssemblyOrExcuse.writeIntervalFile(params.intervalFile, header, intervals, alignedAssemblyOrExcuseList);
        }

        // write alignments of the assembled contigs
        AlignedAssemblyOrExcuse.writeAssemblySAMFile(outputAssemblyAlignments, alignedAssemblyOrExcuseList, header, params.assembliesSortOrder);
        log("Wrote SAM file of aligned contigs.", logger);

        return new AssembledEvidenceResults(evidenceScanResults.readMetadata, intervals, alignedAssemblyOrExcuseList,
                                            evidenceScanResults.evidenceTargetLinks);
    }

    public static final class AssembledEvidenceResults {
        final ReadMetadata readMetadata;
        final List<SVInterval> assembledIntervals;
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList;
        final List<EvidenceTargetLink> evidenceTargetLinks;

        public AssembledEvidenceResults(final ReadMetadata readMetadata,
                                        final List<SVInterval> assembledIntervals,
                                        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList,
                                        final List<EvidenceTargetLink> evidenceTargetLinks) {
            this.readMetadata = readMetadata;
            this.assembledIntervals = assembledIntervals;
            this.alignedAssemblyOrExcuseList = alignedAssemblyOrExcuseList;
            this.evidenceTargetLinks = evidenceTargetLinks;
        }

        public ReadMetadata getReadMetadata() {
            return readMetadata;
        }

        public List<SVInterval> getAssembledIntervals() { return assembledIntervals; }

        public List<AlignedAssemblyOrExcuse> getAlignedAssemblyOrExcuseList() {
            return alignedAssemblyOrExcuseList;
        }

        public List<EvidenceTargetLink> getEvidenceTargetLinks() {
            return evidenceTargetLinks;
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
    private static EvidenceScanResults getMappedQNamesSet(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);
        final List<List<BreakpointEvidence>> externalEvidence =
                readExternalEvidence(params.externalEvidenceFile, readMetadata,
                                        params.externalEvidenceWeight, params.externalEvidenceUncertainty);
        log("External evidence retrieved.", logger);
        final Broadcast<List<List<BreakpointEvidence>>> broadcastExternalEvidence = ctx.broadcast(externalEvidence);
        final Tuple2<List<SVInterval>, List<EvidenceTargetLink>> intervalsAndEvidenceTargetLinks =
                getIntervalsAndEvidenceTargetLinks(params, broadcastMetadata, broadcastExternalEvidence, header,
                        unfilteredReads, filter, logger);
        List<SVInterval> intervals = intervalsAndEvidenceTargetLinks._1();

        SparkUtils.destroyBroadcast(broadcastExternalEvidence, "external evidence");

        final int nIntervals = intervals.size();
        log("Discovered " + nIntervals + " intervals.", logger);

        if ( nIntervals == 0 ) return new EvidenceScanResults(readMetadata, intervals, intervalsAndEvidenceTargetLinks._2(), null);

        if ( params.exclusionIntervalsFile != null ) {
            intervals = removeIntervalsNearGapsAndLog(intervals, params.exclusionIntervalPadding, readMetadata,
                    params.exclusionIntervalsFile, logger);
        }
        intervals = removeHighCoverageIntervalsAndLog(
                                    params, ctx, broadcastMetadata, intervals, unfilteredReads, filter, logger);

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap =
                getQNames(params, ctx, broadcastMetadata, intervals, unfilteredReads, filter);

        SparkUtils.destroyBroadcast(broadcastMetadata, "read metadata");

        if ( params.qNamesMappedFile != null ) {
            QNameAndInterval.writeQNames(params.qNamesMappedFile, qNamesMultiMap);
        }
        log("Discovered " + qNamesMultiMap.size() + " mapped template names.", logger);

        return new EvidenceScanResults(readMetadata, intervals, intervalsAndEvidenceTargetLinks._2(), qNamesMultiMap);
    }

    static final class EvidenceScanResults {
        final ReadMetadata readMetadata;
        final List<SVInterval> intervals;
        final List<EvidenceTargetLink> evidenceTargetLinks;
        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesForAssemblyMultiMap;

        public EvidenceScanResults(final ReadMetadata readMetadata,
                                   final List<SVInterval> intervals,
                                   final List<EvidenceTargetLink> evidenceTargetLinks,
                                   final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesForAssemblyMultiMap) {
            this.readMetadata = readMetadata;
            this.intervals = intervals;
            this.evidenceTargetLinks = evidenceTargetLinks;
            this.qNamesForAssemblyMultiMap = qNamesForAssemblyMultiMap;
        }
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

    @VisibleForTesting
    static List<List<BreakpointEvidence>> readExternalEvidence( final String path,
                                                                final ReadMetadata readMetadata,
                                                                final int externalEvidenceWeight,
                                                                final int externalEvidenceUncertainty ) {
        final int nPartitions = readMetadata.getNPartitions();
        final List<List<BreakpointEvidence>> evidenceByPartition = new ArrayList<>(nPartitions);
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            evidenceByPartition.add(new ArrayList<>());
        }

        if ( path == null ) {
            return evidenceByPartition;
        }

        final SVLocation[] partitionBoundaries = new SVLocation[nPartitions+1];
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            final ReadMetadata.PartitionBounds bounds = readMetadata.getPartitionBounds(idx);
            partitionBoundaries[idx] = new SVLocation(bounds.getFirstContigID(), bounds.getFirstStart());
        }
        partitionBoundaries[nPartitions] = new SVLocation(Integer.MAX_VALUE, Integer.MAX_VALUE);
        for ( int idx = 0; idx != nPartitions; ++idx ) {
            if ( partitionBoundaries[idx].compareTo(partitionBoundaries[idx+1]) > 0 ) {
                throw new GATKException("Partition boundaries are not coordinate sorted.");
            }
        }

        final Map<String, Integer> contigNameMap = readMetadata.getContigNameMap();

        SVLocation prevLocation = new SVLocation(0, 0);
        int partitionIdx = 0;
        try ( final FeatureDataSource<Feature> dataSource = new FeatureDataSource<>(path, null, 0, null) ) {
            for ( final Feature feature : dataSource ) {
                final Integer contigID = contigNameMap.get(feature.getContig());
                if ( contigID == null ) {
                    throw new UserException(path + " contains a contig name not present in the BAM header: " + feature.getContig());
                }
                final int featureStart = feature.getStart();
                final SVLocation featureLocation = new SVLocation(contigID, featureStart);
                if ( prevLocation.compareTo(featureLocation) > 0 ) {
                    throw new UserException("Features in " + path + " are not coordinate sorted.");
                }
                prevLocation = featureLocation;
                while ( featureLocation.compareTo(partitionBoundaries[partitionIdx+1]) >= 0 ) {
                    partitionIdx += 1;
                }
                final List<BreakpointEvidence> partitionEvidence = evidenceByPartition.get(partitionIdx);
                final int featureEnd = feature.getEnd();
                if ( featureEnd - featureStart <= 2*externalEvidenceUncertainty ) {
                    partitionEvidence.add(new ExternalEvidence(contigID, featureStart, featureEnd, externalEvidenceWeight));
                } else {
                    final SVInterval interval1 =
                            BreakpointEvidence.fixedWidthInterval(contigID, featureStart, externalEvidenceUncertainty);
                    partitionEvidence.add(new ExternalEvidence(interval1, externalEvidenceWeight));
                    final SVInterval interval2 =
                            BreakpointEvidence.fixedWidthInterval(contigID, featureEnd, externalEvidenceUncertainty);
                    final SVLocation startLocation2 = interval2.getStartLocation();
                    int partitionIdx2 = partitionIdx;
                    if ( startLocation2.compareTo(partitionBoundaries[partitionIdx+1]) >= 0 ) {
                        partitionIdx2 = Arrays.binarySearch(partitionBoundaries, interval2.getStartLocation());
                        // Binary search returns the complement of the insert position when the exact key is not found.
                        // I.e., ~partitionIdx2 will point to the smallest entry larger than the sought key.
                        // We want the partition before that one, if it exists.
                        if ( partitionIdx2 < 0 ) {
                            partitionIdx2 = ~partitionIdx2;
                            if ( partitionIdx2 > 0 ) partitionIdx2 -= 1;
                        }
                    }
                    evidenceByPartition.get(partitionIdx2).add(new ExternalEvidence(interval2, externalEvidenceWeight));
                }
            }
        }
        return evidenceByPartition;
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
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Tuple2<List<AlignedAssemblyOrExcuse>, HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerAndIntervalsSet(params, readMetadata, ctx, qNamesMultiMap, nIntervals,
                                        unfilteredReads, filter, logger);

        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmersAndIntervals =
                removeUbiquitousKmers(params, readMetadata, ctx, kmerIntervalsAndDispositions._2(), unfilteredReads, filter, logger);

        qNamesMultiMap.addAll(getAssemblyQNames(params, ctx, kmersAndIntervals, unfilteredReads, filter));

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
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger)
    {
        final Set<SVKmer> kmerKillSet =
                SVFileUtils.readKmersFile(params.kmersToIgnoreFile, params.kSize);
        if ( params.adapterSequence != null ) {
            SVKmerizer.stream(params.adapterSequence, params.kSize, 0, new SVKmerLong())
                    .forEach(kmer -> kmerKillSet.add(kmer.canonical(params.kSize)));
        }
        log("Ignoring " + kmerKillSet.size() + " genomically common kmers.", logger);

        final Tuple2<List<AlignedAssemblyOrExcuse>, List<KmerAndInterval>> kmerIntervalsAndDispositions =
                getKmerIntervals(params, readMetadata, ctx, qNamesMultiMap, nIntervals, kmerKillSet,
                                    unfilteredReads, filter, logger);
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

        final int[] counts = new int[nIntervals];
        for ( final QNameAndInterval qNameAndInterval : qNamesMultiMap ) {
            counts[qNameAndInterval.getIntervalId()] += 1;
        }
        final ComplexityPartitioner partitioner = new ComplexityPartitioner(counts);

        final Broadcast<HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> broadcastQNamesMultiMap =
                ctx.broadcast(qNamesMultiMap);
        final List<AlignedAssemblyOrExcuse> intervalDispositions =
            unfilteredReads
                .mapPartitionsToPair(readItr ->
                        new ReadsForQNamesFinder(broadcastQNamesMultiMap.value(), nIntervals,
                                includeMappingLocation, readItr, filter).iterator(), false)
                .combineByKey(x -> x,
                                SVUtils::concatenateLists,
                                SVUtils::concatenateLists,
                                partitioner, false, null)
                .map(localAssemblyHandler::apply)
                .collect();

        SparkUtils.destroyBroadcast(broadcastQNamesMultiMap, "QNames multi map");
        BwaMemIndexCache.closeAllDistributedInstances(ctx);

        return intervalDispositions;
    }

    public static final class IntPair {
        private final int int1;
        private final int int2;

        public IntPair( final int int1, final int int2 ) {
            this.int1 = int1;
            this.int2 = int2;
        }

        public int int1() { return int1; }
        public int int2() { return int2; }

        public static IntPair reduce( final IntPair pair1, final IntPair pair2 ) {
            return new IntPair(pair1.int1 + pair2.int1, pair1.int2 + pair2.int2);
        }
    }

    /**
     * For a set of interesting kmers, count occurrences of each over all reads, and remove those
     * that appear too frequently from the set.
     */
    private static HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> removeUbiquitousKmers(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmersAndIntervals,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger ) {
        final Broadcast<HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> broadcastKmersAndIntervals =
                ctx.broadcast(kmersAndIntervals);

        final int kmersPerPartition = kmersAndIntervals.size();
        final int kSize = params.kSize;
        final int partitionSpan = readMetadata.getMedianPartitionSpan();
        // Especially if reads are long and partitions small, the number of partitions in which we observe a kmer will
        // be larger than the actual number of copies at independent loci.
        // This calculation computes the expected inflator given the typical partition size, read length, and kmer size.
        final float partitionsPerLocus = (partitionSpan + readMetadata.getAvgReadLen() - params.kSize) / partitionSpan;
        final int maxPartitionCount = (int)(params.cleanerMaxCopyNumber * partitionsPerLocus);
        logger.info("Cleanup: maxPartitions=" + maxPartitionCount);
        final int maxCount = (int)(params.cleanerMaxCopyNumber * readMetadata.getAccurateKmerCoverage(kSize));
        // The IntPair class is used to keep two counts:
        //  int1 is the number of partitions in which we observe the kmer.
        //  int2 is the total number of observations of the kmer.
        final List<SVKmer> ubiquitousKmers =
                unfilteredReads
                        .filter(filter::notJunk)
                        .filter(filter::isPrimaryLine)
                        .mapPartitions(readItr ->
                            new KmerCounter(kSize,kmersPerPartition,broadcastKmersAndIntervals.getValue()).apply(readItr))
                        .mapToPair(kmerAndCount -> new Tuple2<>(kmerAndCount.getKey(), new IntPair(1, kmerAndCount.getValue())))
                        .reduceByKey(IntPair::reduce)
                        .filter(pair -> {
                            final IntPair intPair = pair._2();
                            return intPair.int1() > maxPartitionCount || intPair.int2() > maxCount;
                        })
                        .map(Tuple2::_1)
                        .collect();

        for ( final SVKmer kmer : ubiquitousKmers ) {
            final Iterator<KmerAndInterval> entryItr = kmersAndIntervals.findEach(kmer);
            while ( entryItr.hasNext() ) {
                entryItr.next();
                entryItr.remove();
            }
        }

        SparkUtils.destroyBroadcast(broadcastKmersAndIntervals, "kmers and intervals");

        log("Removed "+ubiquitousKmers.size()+" ubiquitous kmers.", logger);

        return kmersAndIntervals;
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
        final Broadcast<HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval>> broadcastKmersAndIntervals =
                ctx.broadcast(kmerMultiMap);

        final int kSize = params.kSize;
        final List<QNameAndInterval> qNamesAndIntervals =
            unfilteredReads
                .filter(filter::notJunk)
                .filter(filter::isPrimaryLine)
                .mapPartitions(readItr ->
                        new FlatMapGluer<>(new QNameIntervalFinder(kSize,broadcastKmersAndIntervals.getValue()), readItr))
                .collect();

        SparkUtils.destroyBroadcast(broadcastKmersAndIntervals, "cleaned kmers and intervals");

        return qNamesAndIntervals;
    }

    /** find kmers for each interval */
    @VisibleForTesting static Tuple2<List<AlignedAssemblyOrExcuse>, List<KmerAndInterval>> getKmerIntervals(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final ReadMetadata readMetadata,
            final JavaSparkContext ctx,
            final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap,
            final int nIntervals,
            final Set<SVKmer> kmerKillSet,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger ) {

        final Broadcast<Set<SVKmer>> broadcastKmerKillSet = ctx.broadcast(kmerKillSet);
        final Broadcast<HopscotchUniqueMultiMap<String, Integer, QNameAndInterval>> broadcastQNameAndIntervalsMultiMap =
                ctx.broadcast(qNamesMultiMap);

        // given a set of template names with interval IDs and a kill set of ubiquitous kmers,
        // produce a set of interesting kmers for each interval ID
        final int kSize = params.kSize;
        final int kmersPerPartition = (int)(2 * readMetadata.getNRefBases() / readMetadata.getNPartitions());
        // diploid coverage from read metadata is:
        //  {number of sequenced bases}/{size of genome}*{a deflation factor that accounts for the number of erroneous
        //  kmers that will be produced given the mean base quality}
        final float effectiveDiploidCoverage = readMetadata.getAccurateKmerCoverage(kSize);
        // haploid coverage is half the diploid coverage
        final float effectiveHaploidCoverage = effectiveDiploidCoverage / 2;
        // expected distribution of the coverage is expected to be poisson with np (==coverage) large enough
        // that it'll look pretty similar to a normal distribution with mean and variance equal to the coverage.
        final float effectiveHaploidCoverageSD = (float)Math.sqrt(effectiveHaploidCoverage);
        // take mean - 3*sigma as the minimum we expect to see by chance
        final float minCoverage = effectiveHaploidCoverage - 3 * effectiveHaploidCoverageSD;
        // but don't go smaller than the specified min count
        final int minKmers = Math.max(params.cleanerMinKmerCount, Math.round(minCoverage));
        // expected number of observations given the copy number and coverage
        final int maxKmers = Math.round(params.cleanerMaxCopyNumber * effectiveDiploidCoverage);
        logger.info("Cleanup: minKmers=" + minKmers + " maxKmers=" + maxKmers);
        final int maxIntervals = params.cleanerMaxIntervals;
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
                        new KmerCleaner(itr, kmersPerPartition, minKmers, maxKmers, maxIntervals).iterator())
                .collect();

        SparkUtils.destroyBroadcast(broadcastQNameAndIntervalsMultiMap, "QNames and intervals");
        SparkUtils.destroyBroadcast(broadcastKmerKillSet, "kmer kill set");

        final int[] intervalKmerCounts = new int[nIntervals];
        for ( final KmerAndInterval kmerAndInterval : kmerIntervals ) {
            intervalKmerCounts[kmerAndInterval.getIntervalId()] += 1;
        }
        final Set<Integer> intervalsToKill = new HashSet<>();
        final List<AlignedAssemblyOrExcuse> intervalDispositions = new ArrayList<>();
        for ( int idx = 0; idx != intervalKmerCounts.length; ++idx ) {
            if ( intervalKmerCounts[idx] < params.minKmersPerInterval ) {
                intervalsToKill.add(idx);
                intervalDispositions.add(new AlignedAssemblyOrExcuse(idx,"FASTQ not written -- too few kmers"));
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
                new TreeSet<>(SVFileUtils.readIntervalsFile(exclusionIntervalsFile, contigNameMap));
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

        SparkUtils.destroyBroadcast(broadcastIntervals, "intervals");

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

        SparkUtils.destroyBroadcast(broadcastIntervals, "intervals");

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNamesMultiMap =
                new HopscotchUniqueMultiMap<>(params.assemblyToMappedSizeRatioGuess*qNameAndIntervalList.size());
        qNamesMultiMap.addAll(qNameAndIntervalList);
        return qNamesMultiMap;
    }

    /**
     * Identify funky reads that support a hypothesis of a breakpoint in the vicinity, group the reads,
     * and declare a breakpoint interval where there is sufficient density of evidence.
     */
    @VisibleForTesting static Tuple2<List<SVInterval>, List<EvidenceTargetLink>> getIntervalsAndEvidenceTargetLinks(
            final FindBreakpointEvidenceSparkArgumentCollection params,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final Broadcast<List<List<BreakpointEvidence>>> broadcastExternalEvidenceByPartition,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final SVReadFilter filter,
            final Logger logger ) {
        // find all breakpoint evidence, then filter for pile-ups
        final int nContigs = header.getSequenceDictionary().getSequences().size();
        final int minEvidenceWeight = params.minEvidenceWeight;
        final int minCoherentEvidenceWeight = params.minCoherentEvidenceWeight;
        final int allowedOverhang = params.allowedShortFragmentOverhang;

        // 1) identify well-mapped reads
        // 2) that look like they support a hypothesis of a breakpoint in the vicinity
        // 3a) filter out those that lack supporting evidence from a sufficient number of other reads, except
        // 3b) pass through everything within a fragment length of partition boundaries
        final JavaRDD<BreakpointEvidence> evidenceRDD = unfilteredReads
                .mapPartitions(readItr -> {
                    final GATKRead sentinel = new SAMRecordToGATKReadAdapter(null);
                    return FlatMapGluer.applyMapFunc(
                            new ReadClassifier(broadcastMetadata.value(), sentinel, allowedOverhang, filter),
                            readItr, sentinel);
                }, true);
        evidenceRDD.cache();

        // record the evidence
        if ( params.unfilteredEvidenceDir != null ) {
            evidenceRDD.map(e -> e.stringRep(broadcastMetadata.getValue(), filter.getMinEvidenceMapQ()))
                       .saveAsTextFile(params.unfilteredEvidenceDir);
        }

        final JavaRDD<EvidenceTargetLink> evidenceTargetLinkJavaRDD = evidenceRDD.mapPartitions(
                itr -> {
                    final ReadMetadata readMetadata = broadcastMetadata.getValue();
                    final EvidenceTargetLinkClusterer clusterer = new EvidenceTargetLinkClusterer(readMetadata,
                            filter.getMinEvidenceMapQ());
                    return clusterer.cluster(itr);
                }).filter(link -> link.readPairs >= 2 || link.splitReads >= 1);

        final List<EvidenceTargetLink> evidenceTargetLinks = evidenceTargetLinkJavaRDD.collect();

        log("Collected " + evidenceTargetLinks.size() + " evidence target links", logger);

        writeTargetLinks(broadcastMetadata, evidenceTargetLinks, params.targetLinkFile);

        final JavaRDD<BreakpointEvidence> filteredEvidenceRDD =
                evidenceRDD
                    .mapPartitionsWithIndex( (idx, evidenceItr1) -> {
                            final ReadMetadata readMetadata = broadcastMetadata.value();
                            final PartitionCrossingChecker xChecker =
                                    new PartitionCrossingChecker(idx, readMetadata,
                                                                 readMetadata.getMaxMedianFragmentSize());
                            final Iterator<BreakpointEvidence> evidenceItr2 =
                                    broadcastExternalEvidenceByPartition.value().get(idx).iterator();
                            final List<Iterator<BreakpointEvidence>> evidenceItrList = new ArrayList<>(2);
                            evidenceItrList.add(evidenceItr1);
                            evidenceItrList.add(evidenceItr2);
                            final Iterator<BreakpointEvidence> evidenceItr =
                                    FlatMapGluer.concatIterators(evidenceItrList.iterator());
                            return new BreakpointDensityFilter(evidenceItr,readMetadata,
                                    minEvidenceWeight,minCoherentEvidenceWeight,xChecker,
                                    filter.getMinEvidenceMapQ());
                        }, true);

        filteredEvidenceRDD.cache();

        // record the evidence
        if ( params.evidenceDir != null ) {
            filteredEvidenceRDD
                    .filter(BreakpointEvidence::isValidated)
                    .saveAsTextFile(params.evidenceDir);
        }

        // replace clumps of evidence with a single aggregate indicator of evidence (to save memory), except
        // don't clump anything within two fragment lengths of a partition boundary.
        // collect the whole mess in the driver.
        final int maxFragmentSize = broadcastMetadata.value().getMaxMedianFragmentSize();
        final List<BreakpointEvidence> collectedEvidence =
                filteredEvidenceRDD
                        .mapPartitionsWithIndex( (idx, readEvidenceItr) ->
                                new FlatMapGluer<>(
                                        new BreakpointEvidenceClusterer(maxFragmentSize,
                                                new PartitionCrossingChecker(idx,broadcastMetadata.value(),2*maxFragmentSize)),
                                        readEvidenceItr,
                                        new BreakpointEvidence(new SVInterval(nContigs,1,1),0,false)), true)
                        .collect();

        filteredEvidenceRDD.unpersist();
        evidenceRDD.unpersist();

        // reapply the density filter (all data collected -- no more worry about partition boundaries).
        final Iterator<BreakpointEvidence> evidenceIterator =
                new BreakpointDensityFilter(collectedEvidence.iterator(),
                        broadcastMetadata.value(), minEvidenceWeight, minCoherentEvidenceWeight,
                        new PartitionCrossingChecker(), filter.getMinEvidenceMapQ());
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
        final Iterator<BreakpointEvidence> evidenceIterator2 =
            new FlatMapGluer<>(new BreakpointEvidenceClusterer(maxFragmentSize,new PartitionCrossingChecker()),
                               allEvidence.iterator(),
                               new BreakpointEvidence(new SVInterval(nContigs,1,1),0,false));

        final List<SVInterval> intervals = new ArrayList<>(allEvidence.size());
        while ( evidenceIterator2.hasNext() ) {
            intervals.add(evidenceIterator2.next().getLocation());
        }

        return new Tuple2<>(intervals, evidenceTargetLinks);
    }

    private static void writeTargetLinks(final Broadcast<ReadMetadata> broadcastMetadata,
                                         final List<EvidenceTargetLink> targetLinks, final String targetLinkFile) {
        if ( targetLinkFile != null ) {
            try (final OutputStreamWriter writer =
                         new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(targetLinkFile)))) {
                targetLinks.iterator().forEachRemaining(entry -> {
                    final String bedpeRecord = entry.toBedpeString(broadcastMetadata.getValue());
                    try {
                        writer.write(bedpeRecord + "\n");
                    } catch (final IOException ioe) {
                        throw new GATKException("Can't write target links to "+ targetLinkFile, ioe);
                    }
                });
            } catch ( final IOException ioe ) {
                throw new GATKException("Can't write target links to "+ targetLinkFile, ioe);
            }
        }
    }

    private static void log(final String message, final Logger logger) {
        logger.info(message);
    }
}
