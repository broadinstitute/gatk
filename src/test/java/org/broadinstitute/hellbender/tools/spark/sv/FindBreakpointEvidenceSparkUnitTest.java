package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchUniqueMultiMap;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchSet;
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

    private final FindBreakpointEvidenceSpark.Params params = FindBreakpointEvidenceSpark.defaultParams;
    private final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
    private final ReadsSparkSource readsSource = new ReadsSparkSource(ctx);
    private final SAMFileHeader header = readsSource.getHeader(readsFile, null);
    private final JavaRDD<GATKRead> reads = readsSource.getParallelReads(readsFile, null, null, 0L);
    private final JavaRDD<GATKRead> mappedReads = reads.filter(read -> !read.isUnmapped());
    private final ReadMetadata readMetadataExpected = new ReadMetadata(header, reads);
    private final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadataExpected);
    private final FindBreakpointEvidenceSpark.Locations locations =
        new FindBreakpointEvidenceSpark.Locations(null, null, null, null, null, null, null);
    private final Set<String> expectedQNames = loadExpectedQNames(qNamesFile);
    private final Set<String> expectedAssemblyQNames = loadExpectedQNames(asmQNamesFile);
    private final List<SVInterval> expectedIntervalList = Arrays.asList(testIntervals);

    @Test(groups = "spark")
    public void getIntervalsTest() {
        final List<SVInterval> actualIntervals =
                FindBreakpointEvidenceSpark.getIntervals(params, broadcastMetadata, header, mappedReads, locations);
        Assert.assertEquals(actualIntervals, expectedIntervalList);
    }

    @Test(groups = "spark")
    public void getQNamesTest() {
        final Set<String> actualQNames = new HashSet<>();
        FindBreakpointEvidenceSpark.getQNames(params, ctx, broadcastMetadata, expectedIntervalList, mappedReads)
                .stream()
                .map(FindBreakpointEvidenceSpark.QNameAndInterval::getKey)
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

        final HopscotchUniqueMultiMap<String, Integer, FindBreakpointEvidenceSpark.QNameAndInterval> qNameMultiMap =
                new HopscotchUniqueMultiMap<>(expectedAssemblyQNames.size());

        // an empty qname map should produce a "too few kmers" disposition for the interval
        final Map<Integer, String> dispositions =
                FindBreakpointEvidenceSpark.getKmerIntervals(params, ctx, qNameMultiMap, 1, Collections.emptySet(), reads, locations, null)._1();
        Assert.assertEquals(dispositions.size(), 1);
        Assert.assertTrue(dispositions.get(0).contains("too few"));

        expectedQNames.stream()
                .map(qName -> new FindBreakpointEvidenceSpark.QNameAndInterval(qName, 0))
                .forEach(qNameMultiMap::add);
        final HopscotchUniqueMultiMap<SVKmer, Integer, FindBreakpointEvidenceSpark.KmerAndInterval> actualKmerAndIntervalSet =
                new HopscotchUniqueMultiMap<>(
                        FindBreakpointEvidenceSpark.getKmerIntervals(params, ctx, qNameMultiMap, 1, new HopscotchSet<>(0),
                                reads, locations, null)._2());
        final Set<SVKmer> expectedKmers = SVUtils.readKmersFile(params.kSize, kmersFile, null, kmer);
        Assert.assertEquals(actualKmerAndIntervalSet.size(), expectedKmers.size());
        for ( final FindBreakpointEvidenceSpark.KmerAndInterval kmerAndInterval : actualKmerAndIntervalSet ) {
            Assert.assertTrue(expectedKmers.contains(kmerAndInterval.getKey()));
        }
    }

    @Test(groups = "spark")
    public void getAssemblyQNamesTest() throws FileNotFoundException {
        final Set<SVKmer> expectedKmers = SVUtils.readKmersFile(params.kSize, kmersFile, null, new SVKmerLong(params.kSize));
        final HopscotchUniqueMultiMap<SVKmer, Integer, FindBreakpointEvidenceSpark.KmerAndInterval> kmerAndIntervalSet =
                new HopscotchUniqueMultiMap<>(expectedKmers.size());
        expectedKmers.stream().
                map(kmer -> new FindBreakpointEvidenceSpark.KmerAndInterval(kmer, 0))
                .forEach(kmerAndIntervalSet::add);
        final Set<String> actualAssemblyQNames = new HashSet<>();
        FindBreakpointEvidenceSpark.getAssemblyQNames(params, ctx, kmerAndIntervalSet, reads)
                .stream()
                .map(FindBreakpointEvidenceSpark.QNameAndInterval::getKey)
                .forEach(actualAssemblyQNames::add);
        Assert.assertEquals(actualAssemblyQNames, expectedAssemblyQNames);
    }

    @Test(groups = "spark")
    public void generateFastqsTest() {
        final HopscotchUniqueMultiMap<String, Integer, FindBreakpointEvidenceSpark.QNameAndInterval> qNameMultiMap =
                new HopscotchUniqueMultiMap<>(expectedAssemblyQNames.size());
        expectedAssemblyQNames.stream()
                .map(qName -> new FindBreakpointEvidenceSpark.QNameAndInterval(qName, 0))
                .forEach(qNameMultiMap::add);
        final String expectedFile = fastqFile; // this is NOT redundant -- local copy necessary for lambda serialization
        FindBreakpointEvidenceSpark.generateFastqs(ctx, qNameMultiMap, reads, 2, true,
                intervalAndFastqBytes -> compareFastqs(intervalAndFastqBytes, expectedFile));
    }

    private static Tuple2<Integer, String> compareFastqs(
            final Tuple2<Integer, List<byte[]>> intervalAndFastqBytes,
            final String fastqFile ) {
        final List<byte[]> fastqList = intervalAndFastqBytes._2;
        SVFastqUtils.sortFastqRecords(fastqList);
        final byte[] concatenatedFastqs = new byte[fastqList.stream().mapToInt(fastq -> fastq.length).sum()];
        int idx = 0;
        for ( final byte[] fastq : fastqList ) {
            System.arraycopy(fastq, 0, concatenatedFastqs, idx, fastq.length);
            idx += fastq.length;
        }

        final String expectedFile = fastqFile+intervalAndFastqBytes._1();
        final ByteArrayInputStream actualStream = new ByteArrayInputStream(concatenatedFastqs);
        try( InputStream expectedStream = new BufferedInputStream(new FileInputStream(expectedFile)) ) {
            int val;
            while ( (val = expectedStream.read()) != -1 )
                Assert.assertEquals(actualStream.read(), val);
            Assert.assertEquals(actualStream.read(), -1);
        }
        catch ( final IOException ioe ) {
            throw new GATKException("can't read expected values", ioe);
        }
        return new Tuple2<>(intervalAndFastqBytes._1(), "hello");
    }

    private Set<String> loadExpectedQNames( final String fileName ) {
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
