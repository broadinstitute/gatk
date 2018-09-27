package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.IntHistogramTest;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;

public final class FindBreakpointEvidenceSparkUnitTest extends GATKBaseTest {
    private static final SVInterval[] testIntervals =
            { new SVInterval(1, 43349482, 43350671), new SVInterval(1, 43353045, 43353870) };

    private final String readsFile = largeFileTestDir + "/sv/SVIntegrationTest_hg19.bam";
    private final String toolDir = getToolTestDataDir();
    private final String qNamesFile = toolDir+"SVBreakpointsTest.qnames";
    private final String kmersFile = toolDir+"SVBreakpointsTest.kmers";
    private final String asmQNamesFile = toolDir+"SVBreakpointsTest.assembly.qnames";
    private final String fastqFile = toolDir+"SVBreakpointsTest.assembly.";

    private final FindBreakpointEvidenceSparkArgumentCollection params =
            new FindBreakpointEvidenceSparkArgumentCollection();
    private final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
    private final ReadsSparkSource readsSource = new ReadsSparkSource(ctx);
    private final SAMFileHeader header = readsSource.getHeader(readsFile, null);
    private final JavaRDD<GATKRead> reads = readsSource.getParallelReads(readsFile, null, null, 0L);
    private final SVReadFilter filter = new SVReadFilter(params);
    private final ReadMetadata readMetadataExpected =
            new ReadMetadata(Collections.emptySet(), header,
                                new LibraryStatistics(new IntHistogram.CDF(IntHistogramTest.genLogNormalSample(320, 129, 10000)),
                                        60000000000L, 600000000L, 1200000000000L, 3000000000L),
                new ReadMetadata.PartitionBounds[]
                        { new ReadMetadata.PartitionBounds(0, 1, 1, 10000, 9999)},
                    100, 10, 30);
    private final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadataExpected);
    private final Broadcast<SVIntervalTree<SVInterval>> broadcastRegionsToIgnore = ctx.broadcast(new SVIntervalTree<>());
    private final List<List<BreakpointEvidence>> externalEvidence =
            FindBreakpointEvidenceSpark.readExternalEvidence(null, readMetadataExpected,
                                                    params.externalEvidenceWeight, params.externalEvidenceUncertainty);
    private final Broadcast<List<List<BreakpointEvidence>>> broadcastExternalEvidence = ctx.broadcast(externalEvidence);
    private final Set<String> expectedQNames = loadExpectedQNames(qNamesFile);
    private final Set<String> expectedAssemblyQNames = loadExpectedQNames(asmQNamesFile);
    private final List<SVInterval> expectedIntervalList = Arrays.asList(testIntervals);

    @Test(groups = "sv")
    public void getIntervalsTest() {
        // Now that thresholds scale with coverage, changing coverage alters thresholds. Re-fix thresholds to match
        // original test:
        final FindBreakpointEvidenceSparkArgumentCollection intervalsParams =
                new FindBreakpointEvidenceSparkArgumentCollection();
        intervalsParams.minEvidenceWeightPerCoverage = 15.0 / broadcastMetadata.getValue().getCoverage();
        intervalsParams.minCoherentEvidenceWeightPerCoverage = 7.0 / broadcastMetadata.getValue().getCoverage();
        // run test as previously:
        final List<SVInterval> actualIntervals =
                FindBreakpointEvidenceSpark.getIntervalsAndEvidenceTargetLinks(intervalsParams,broadcastMetadata,
                        broadcastExternalEvidence,header,reads,filter,logger, broadcastRegionsToIgnore)._1();
        Assert.assertEquals(actualIntervals, expectedIntervalList);
    }

    @Test(groups = "sv")
    public void getQNamesTest() {
        final Set<String> actualQNames = new HashSet<>();
        FindBreakpointEvidenceSpark.getQNames(params, ctx, broadcastMetadata, expectedIntervalList, reads, filter, broadcastRegionsToIgnore)
                .stream()
                .map(QNameAndInterval::getKey)
                .forEach(actualQNames::add);
        Assert.assertEquals(actualQNames, expectedQNames);
    }

    @Test(groups = "sv")
    public void getKmerIntervalsTest() {
        final SVKmer kmer = new SVKmerLong(params.kSize);
        final Set<SVKmer> killSet = new HashSet<>();
        final String seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG";
        killSet.add(SVKmerizer.toKmer(seq1,kmer));
        final String seq2 = "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
        killSet.add(SVKmerizer.toKmer(seq2,kmer));
        Assert.assertEquals( killSet.size(), 2);

        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNameMultiMap =
                new HopscotchUniqueMultiMap<>(expectedAssemblyQNames.size());

        // an empty qname map should produce a "too few kmers" disposition for the interval
        final List<AlignedAssemblyOrExcuse> alignedAssemblyOrExcuseList =
                FindBreakpointEvidenceSpark.getKmerIntervals(
                        params, readMetadataExpected, ctx, qNameMultiMap, 1,
                        Collections.emptySet(), reads, filter, logger)._1();
        Assert.assertEquals(alignedAssemblyOrExcuseList.size(), 1);
        Assert.assertTrue(alignedAssemblyOrExcuseList.get(0).getErrorMessage().contains("too few"));

        expectedQNames.stream()
                .map(qName -> new QNameAndInterval(qName, 0))
                .forEach(qNameMultiMap::add);
        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> actualKmerAndIntervalSet =
                new HopscotchUniqueMultiMap<>(
                        FindBreakpointEvidenceSpark.getKmerIntervals(params, readMetadataExpected, ctx, qNameMultiMap,
                                1, new HopscotchSet<>(0),
                                reads, filter, logger)._2());
        final Set<SVKmer> actualKmers = new HashSet<>(SVUtils.hashMapCapacity(actualKmerAndIntervalSet.size()));
        for ( final KmerAndInterval kmerAndInterval : actualKmerAndIntervalSet ) {
            actualKmers.add(kmerAndInterval.getKey());
        }
        final Set<SVKmer> expectedKmers = SVFileUtils.readKmersFile(kmersFile, params.kSize);
        Assert.assertEquals(actualKmers, expectedKmers);
    }

    @Test(groups = "sv")
    public void getAssemblyQNamesTest() {
        final Set<SVKmer> expectedKmers = SVFileUtils.readKmersFile(kmersFile, params.kSize);
        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> kmerAndIntervalSet =
                new HopscotchUniqueMultiMap<>(expectedKmers.size());
        expectedKmers.stream().
                map(kmer -> new KmerAndInterval(kmer, 0))
                .forEach(kmerAndIntervalSet::add);
        final Set<String> actualAssemblyQNames = new HashSet<>();
        FindBreakpointEvidenceSpark.getAssemblyQNames(params, ctx, kmerAndIntervalSet, reads, filter)
                .stream()
                .map(QNameAndInterval::getKey)
                .forEach(actualAssemblyQNames::add);
        Assert.assertEquals(actualAssemblyQNames, expectedAssemblyQNames);
    }

    @Test(groups = "sv")
    public void generateFastqsTest() {
        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNameMultiMap =
                new HopscotchUniqueMultiMap<>(expectedAssemblyQNames.size());
        expectedAssemblyQNames.stream()
                .map(qName -> new QNameAndInterval(qName, 0))
                .forEach(qNameMultiMap::add);
        FindBreakpointEvidenceSpark.handleAssemblies(ctx,qNameMultiMap,reads,filter,2,true,new LocalAssemblyComparator(fastqFile));
    }

    @Test(groups = "sv")
    public void readExternalEvidenceTest() {
        final int evidenceWeight = params.externalEvidenceWeight;
        final int evidenceSlop = params.externalEvidenceUncertainty;
        final ReadMetadata.PartitionBounds bounds = readMetadataExpected.getPartitionBounds(0);
        final SVInterval interval1 = new SVInterval(bounds.getFirstContigID(), bounds.getFirstStart(), bounds.getFirstStart()+100);
        final SVInterval interval2 = new SVInterval(bounds.getLastContigID(), bounds.getLastEnd(), bounds.getLastEnd()+200);
        final File file = createTempFile("test", ".bed");
        try ( final FileWriter writer = new FileWriter(file) ) {
            writer.write(readMetadataExpected.getContigName(interval1.getContig()) + "\t" + (interval1.getStart()-1) + "\t" + interval1.getEnd() + "\n");
            writer.write(readMetadataExpected.getContigName(interval2.getContig()) + "\t" + (interval2.getStart()-1) + "\t" + interval2.getEnd() + "\n");
        } catch ( final IOException ioe ) {
            throw new GATKException("failed to write test bed file", ioe);
        }
        final List<List<BreakpointEvidence>> actual;
        try {
            actual = FindBreakpointEvidenceSpark.readExternalEvidence(file.getPath(), readMetadataExpected,
                                                                        evidenceWeight, evidenceSlop);
        } finally {
            file.delete();
        }
        Assert.assertEquals(actual.size(), readMetadataExpected.getNPartitions());
        final List<BreakpointEvidence> partition0expected = new ArrayList<>(2);
        partition0expected.add(new BreakpointEvidence.ExternalEvidence(interval1, evidenceWeight));
        partition0expected.add(new BreakpointEvidence.ExternalEvidence(interval2, evidenceWeight));
        final List<BreakpointEvidence> partition0actual = actual.get(0);
        Assert.assertEquals(partition0actual.size(), partition0expected.size());
        for ( int idx = 0; idx != partition0actual.size(); ++idx ) {
            Assert.assertTrue(partition0actual.get(idx).equalFields(partition0expected.get(idx)));
        }
    }

    /** This LocalAssemblyHandler compares an assembly with expected results. */
    private static final class LocalAssemblyComparator implements FindBreakpointEvidenceSpark.LocalAssemblyHandler {
        private static final long serialVersionUID = 1L;
        private final String expectedFastqFilePrefix;

        public LocalAssemblyComparator( final String expectedFastqFilePrefix ) { this.expectedFastqFilePrefix = expectedFastqFilePrefix; }

        @Override
        public AlignedAssemblyOrExcuse
                    apply( final Tuple2<Integer,List<SVFastqUtils.FastqRead>> intervalAndFastqBytes ) {
            final List<SVFastqUtils.FastqRead> fastqList = intervalAndFastqBytes._2();
            fastqList.sort(Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
            final File outputFile = createTempFile("output", ".fastq");
            SVFastqUtils.writeFastqFile(outputFile.getAbsolutePath() , fastqList.iterator());

            try ( final BufferedReader actualReader = new BufferedReader(new FileReader(outputFile));
                  final BufferedReader expectedReader =
                          new BufferedReader(new FileReader(expectedFastqFilePrefix + intervalAndFastqBytes._1())) ) {
                final List<String> actualLines = actualReader.lines().collect(Collectors.toList());
                final List<String> expectedLines = expectedReader.lines().collect(Collectors.toList());
                Assert.assertEquals(actualLines.size(), expectedLines.size(), "different number of lines");
                for (int i = 0; i < actualLines.size(); i++) {
                    Assert.assertEquals(actualLines.get(i), expectedLines.get(i), "difference in fastq line " + (i + 1));
                }
            } catch (final IOException ex) {
                throw new GATKException("test problems reading FASTQ's", ex);
            } finally {
                try { outputFile.delete(); } catch (final Exception e) {};
            }
            return new AlignedAssemblyOrExcuse(intervalAndFastqBytes._1(), "hello");
        }
    }

    private static Set<String> loadExpectedQNames( final String fileName ) {
        final HashSet<String> qNameSet = new HashSet<>();
        try( BufferedReader rdr = new BufferedReader(new FileReader(fileName)) ) {
            for ( String line = rdr.readLine(); line != null; line = rdr.readLine() ) {
                qNameSet.add(line);
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Unable to read " + fileName);
        }
        return qNameSet;
    }
}
