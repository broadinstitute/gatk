package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
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

public final class FindBreakpointEvidenceSparkUnitTest extends BaseTest {
    private static final SVInterval[] testIntervals =
            { new SVInterval(1, 33140934, 33141485), new SVInterval(1, 33143116, 33143539) };

    private final String toolDir = getToolTestDataDir();
    private final String readsFile = toolDir+"SVBreakpointsTest.bam";
    private final String qNamesFile = toolDir+"SVBreakpointsTest.qnames";
    private final String kmersFile = toolDir+"SVBreakpointsTest.kmers";
    private final String asmQNamesFile = toolDir+"SVBreakpointsTest.assembly.qnames";
    private final String fastqFile = toolDir+"SVBreakpointsTest.assembly.";

    private final FindBreakpointEvidenceSpark.Params params =
            StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection.defaultParams;
    private final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
    private final ReadsSparkSource readsSource = new ReadsSparkSource(ctx);
    private final SAMFileHeader header = readsSource.getHeader(readsFile, null);
    private final JavaRDD<GATKRead> reads = readsSource.getParallelReads(readsFile, null, null, 0L);
    private final JavaRDD<GATKRead> mappedReads = reads.filter(read -> !read.isUnmapped());
    private final ReadMetadata readMetadataExpected = new ReadMetadata(Collections.emptySet(), header, reads);
    private final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadataExpected);
    private final FindBreakpointEvidenceSpark.Locations locations =
        new FindBreakpointEvidenceSpark.Locations(null, null, null, null, null, null, null, null, null);
    private final Set<String> expectedQNames = loadExpectedQNames(qNamesFile);
    private final Set<String> expectedAssemblyQNames = loadExpectedQNames(asmQNamesFile);
    private final List<SVInterval> expectedIntervalList = Arrays.asList(testIntervals);

    @Test(groups = "spark")
    public void getIntervalsTest() {
        final List<SVInterval> actualIntervals =
                FindBreakpointEvidenceSpark.getIntervals(params,broadcastMetadata,header,mappedReads,locations);
        Assert.assertEquals(actualIntervals, expectedIntervalList);
    }

    @Test(groups = "spark")
    public void getQNamesTest() {
        final Set<String> actualQNames = new HashSet<>();
        FindBreakpointEvidenceSpark.getQNames(params, ctx, broadcastMetadata, expectedIntervalList, mappedReads)
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
                FindBreakpointEvidenceSpark.getKmerIntervals(params, ctx, qNameMultiMap, 1, Collections.emptySet(), reads, locations)._1();
        Assert.assertEquals(alignedAssemblyOrExcuseList.size(), 1);
        Assert.assertTrue(alignedAssemblyOrExcuseList.get(0).getErrorMessage().contains("too few"));

        expectedQNames.stream()
                .map(qName -> new QNameAndInterval(qName, 0))
                .forEach(qNameMultiMap::add);
        final HopscotchUniqueMultiMap<SVKmer, Integer, KmerAndInterval> actualKmerAndIntervalSet =
                new HopscotchUniqueMultiMap<>(
                        FindBreakpointEvidenceSpark.getKmerIntervals(params, ctx, qNameMultiMap, 1, new HopscotchSet<>(0),
                                reads, locations)._2());
        final Set<SVKmer> expectedKmers = SVUtils.readKmersFile(params.kSize, kmersFile, kmer);
        Assert.assertEquals(actualKmerAndIntervalSet.size(), expectedKmers.size());
        for ( final KmerAndInterval kmerAndInterval : actualKmerAndIntervalSet ) {
            Assert.assertTrue(expectedKmers.contains(kmerAndInterval.getKey()));
        }
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
        FindBreakpointEvidenceSpark.getAssemblyQNames(params, ctx, kmerAndIntervalSet, reads)
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
        FindBreakpointEvidenceSpark.handleAssemblies(ctx,qNameMultiMap,reads,2,true,true,new LocalAssemblyComparator(fastqFile));
    }

    /** This LocalAssemblyHandler compares an assembly with expected results. */
    private static final class LocalAssemblyComparator implements FindBreakpointEvidenceSpark.LocalAssemblyHandler {
        private static final long serialVersionUID = 1L;
        private final String expectedFastqFile;

        public LocalAssemblyComparator( final String expectedFastqFile ) { this.expectedFastqFile = expectedFastqFile; }

        @Override
        public AlignedAssemblyOrExcuse
                    apply( final Tuple2<Integer,List<SVFastqUtils.FastqRead>> intervalAndFastqBytes ) {
            final List<SVFastqUtils.FastqRead> fastqList = intervalAndFastqBytes._2();
            fastqList.sort(Comparator.comparing(SVFastqUtils.FastqRead::getName));
            final ByteArrayOutputStream os = new ByteArrayOutputStream();
            try { SVFastqUtils.writeFastqStream(os, fastqList.iterator()); }
            catch ( final IOException ioe ) { throw new GATKException("can't stream fastqs into memory", ioe); }
            final byte[] concatenatedFastqs = os.toByteArray();
            final String expectedFile = expectedFastqFile+intervalAndFastqBytes._1();
            final ByteArrayInputStream actualStream = new ByteArrayInputStream(concatenatedFastqs);
            try( InputStream expectedStream = new BufferedInputStream(new FileInputStream(expectedFile)) ) {
                int val;
                while ( (val = expectedStream.read()) != -1 )
                    Assert.assertEquals(actualStream.read(), val);
                Assert.assertEquals(actualStream.read(), -1);
            }
            catch ( final IOException ioe2 ) {
                throw new GATKException("can't read expected values", ioe2);
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
