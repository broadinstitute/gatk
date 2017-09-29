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
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;

public final class FindBreakpointEvidenceSparkUnitTest extends BaseTest {
    private static final SVInterval[] testIntervals =
            { new SVInterval(1, 33140717, 33141485), new SVInterval(1, 33143109, 33143539) };

    private final String toolDir = getToolTestDataDir();
    private final String readsFile = toolDir+"SVBreakpointsTest.bam";
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
            new ReadMetadata(Collections.emptySet(), header, params.maxTrackedFragmentLength, reads, filter);
    private final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadataExpected);
    private final Set<String> expectedQNames = loadExpectedQNames(qNamesFile);
    private final Set<String> expectedAssemblyQNames = loadExpectedQNames(asmQNamesFile);
    private final List<SVInterval> expectedIntervalList = Arrays.asList(testIntervals);

    @Test(groups = "spark")
    public void getIntervalsTest() {
        final List<SVInterval> actualIntervals =
                FindBreakpointEvidenceSpark.getIntervals(params,broadcastMetadata,header,reads,filter);
        Assert.assertEquals(actualIntervals, expectedIntervalList);
    }

    @Test(groups = "spark")
    public void getQNamesTest() {
        final Set<String> actualQNames = new HashSet<>();
        FindBreakpointEvidenceSpark.getQNames(params, ctx, broadcastMetadata, expectedIntervalList, reads, filter)
                .stream()
                .map(QNameAndInterval::getKey)
                .forEach(actualQNames::add);
        Assert.assertEquals(actualQNames, expectedQNames);
    }

    @Test(groups = "spark")
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
                        params, ctx, qNameMultiMap, 1, Collections.emptySet(), reads, filter)._1();
        Assert.assertEquals(alignedAssemblyOrExcuseList.size(), 1);
        Assert.assertTrue(alignedAssemblyOrExcuseList.get(0).getErrorMessage().contains("too few"));

        expectedQNames.stream()
                .map(qName -> new QNameAndInterval(qName, 0))
                .forEach(qNameMultiMap::add);
        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> actualKmerAndIntervalSet =
                new HopscotchUniqueMultiMap<>(
                        FindBreakpointEvidenceSpark.getKmerIntervals(params, ctx, qNameMultiMap, 1, new HopscotchSet<>(0),
                                reads, filter)._2());
        final Set<SVKmer> actualKmers = new HashSet<SVKmer>(SVUtils.hashMapCapacity(actualKmerAndIntervalSet.size()));
        for ( final KmerAndInterval kmerAndInterval : actualKmerAndIntervalSet ) {
            actualKmers.add(kmerAndInterval.getKey());
        }
        final Set<SVKmer> expectedKmers = SVUtils.readKmersFile(params.kSize, kmersFile, kmer);
        Assert.assertEquals(actualKmers, expectedKmers);
    }

    @Test(groups = "spark")
    public void getAssemblyQNamesTest() throws FileNotFoundException {
        final Set<SVKmer> expectedKmers = SVUtils.readKmersFile(params.kSize, kmersFile, new SVKmerLong(params.kSize));
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

    @Test(groups = "spark")
    public void generateFastqsTest() {
        final HopscotchUniqueMultiMap<String, Integer, QNameAndInterval> qNameMultiMap =
                new HopscotchUniqueMultiMap<>(expectedAssemblyQNames.size());
        expectedAssemblyQNames.stream()
                .map(qName -> new QNameAndInterval(qName, 0))
                .forEach(qNameMultiMap::add);
        FindBreakpointEvidenceSpark.handleAssemblies(ctx,qNameMultiMap,reads,filter,2,true,new LocalAssemblyComparator(fastqFile));
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
            Collections.sort(fastqList, Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
            final File outputFile = createTempFile("output", ".fastq");
            try {
                SVFastqUtils.writeFastqFile(outputFile.getAbsolutePath() , fastqList.iterator());
            } catch (final RuntimeException ex) {
                throw ex;
            }
            try (final BufferedReader actualReader = new BufferedReader(new FileReader(outputFile));
                 final BufferedReader expectedReader = new BufferedReader(new FileReader(expectedFastqFilePrefix + intervalAndFastqBytes._1()))){
                final List<String> actualLines = actualReader.lines().collect(Collectors.toList());
                final List<String> expectedLines = expectedReader.lines().collect(Collectors.toList());
                Assert.assertEquals(actualLines.size(), expectedLines.size(), "different number of lines");
                for (int i = 0; i < actualLines.size(); i++) {
                    Assert.assertEquals(actualLines.get(i), expectedLines.get(i), "difference in fastq line " + (i + 1));
                }
            } catch (final IOException ex) {
                throw new GATKException("test problems", ex);
            } finally {
                try { outputFile.delete(); } catch (final Throwable t) {};
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
