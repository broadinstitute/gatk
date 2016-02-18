package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * Unit tests for FindSVBreakpointsSpark.
 */
public class FindSVBreakpointsSparkUnitTest extends BaseTest {

    private static final String REFERENCE_FILE_NAME = hg19MiniReference;
    private static final List<GATKRead> TEST_READS;
    static {
        TEST_READS = new ArrayList<>(4);
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        TEST_READS.addAll(ArtificialReadUtils.createPair(samHeader, "firstReadPair", 101, 1010, 1382, false, false));
        TEST_READS.addAll(ArtificialReadUtils.createPair(samHeader, "secondReadPair", 101, 1210, 1582, false, false));
        for ( int idx = 0; idx != 4; ++idx ) {
            final GATKRead read = TEST_READS.get(idx);
            read.getBases()[0] = BaseUtils.BASES[idx];
            read.getBaseQualities()[0] = (byte)(30 + idx);
            read.setMappingQuality(FindSVBreakpointsSpark.BreakpointClusterer.MIN_MAPQ);
        }
        TEST_READS.get(0).setCigar("51S50M");
        TEST_READS.get(2).setIsProperlyPaired(false);
        TEST_READS.get(2).setMateIsReverseStrand(false);
        TEST_READS.get(3).setIsProperlyPaired(false);
        TEST_READS.get(3).setIsReverseStrand(false);
    }
    @Test(groups = "spark")
    void kmersToIgnoreTest() throws IOException {
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        final ReferenceMultiSource ref = new ReferenceMultiSource((PipelineOptions)null,
                REFERENCE_FILE_NAME,
                ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final List<SVKmer> badKmers = FindBadGenomicKmersSpark.findBadGenomicKmers(ctx, ref, null, null);
        final File kmersFile = File.createTempFile("kmerKillList", "txt");
        kmersFile.deleteOnExit();
        FindBadGenomicKmersSpark.writeKmersToOutput(new FileOutputStream(kmersFile), badKmers);
        final Set<SVKmer> kmerSet = FindSVBreakpointsSpark.readKmersToIgnore(kmersFile);
        if ( !kmersFile.delete() )
            throw new GATKException("Unable to delete file "+kmersFile);
        Assert.assertEquals(badKmers.size(), kmerSet.size());
        for ( final SVKmer kmer : badKmers )
            Assert.assertTrue(kmerSet.contains(kmer));
    }

    @Test(groups = "spark")
    void writeFastqTest() {
        final File fastqName = FindSVBreakpointsSpark.writeFastq(0L, TEST_READS.iterator(), ".");
        fastqName.deleteOnExit();

        final FastqReader reader = new FastqReader(fastqName);
        final List<FastqRecord> fastqRecords = new ArrayList<>(4);
        while ( reader.hasNext() ) {
            fastqRecords.add(reader.next());
        }

        if ( !fastqName.delete() )
            throw new GATKException("Unable to delete "+fastqName);

        Assert.assertEquals(TEST_READS.size(), fastqRecords.size());
        for ( int idx = 0; idx != TEST_READS.size(); ++idx ) {
            final GATKRead read = TEST_READS.get(idx);
            final FastqRecord fastqRecord = fastqRecords.get(idx);
            Assert.assertEquals(read.getBasesString(), fastqRecord.getReadString());
            Assert.assertEquals(read.getBaseQualities(), SAMUtils.fastqToPhred(fastqRecord.getBaseQualityString()));
        }
    }

    @Test(groups = "spark")
    void eventLocusTest() {
        final FindSVBreakpointsSpark.EventLocus eventLocus1 = new FindSVBreakpointsSpark.EventLocus(100, 3, 0L);
        Assert.assertEquals(eventLocus1.getLocusStart(), 100);
        Assert.assertEquals(eventLocus1.getLocusEnd(), 103);
        final FindSVBreakpointsSpark.EventLocus eventLocus2 = new FindSVBreakpointsSpark.EventLocus(100, 3, 0L);
        Assert.assertEquals(eventLocus1.compareTo(eventLocus2)==0, true);
        final FindSVBreakpointsSpark.EventLocus eventLocus3 = new FindSVBreakpointsSpark.EventLocus(100, 3, 1L);
        Assert.assertEquals(eventLocus1.compareTo(eventLocus3)<0, true);
        final FindSVBreakpointsSpark.EventLocus eventLocus4 = new FindSVBreakpointsSpark.EventLocus(99, 4, 0L);
        Assert.assertEquals(eventLocus1.compareTo(eventLocus4)>0, true);
    }

    @Test(groups = "spark")
    void testSplitReadDetector() {
        final FindSVBreakpointsSpark.SplitReadDetector detector = new FindSVBreakpointsSpark.SplitReadDetector();
        Assert.assertNotNull(detector.getLocusIfReadIsSplit(TEST_READS.get(0),0L));
        Assert.assertNull(detector.getLocusIfReadIsSplit(TEST_READS.get(1),0L));
        Assert.assertNull(detector.getLocusIfReadIsSplit(TEST_READS.get(2),0L));
        Assert.assertNull(detector.getLocusIfReadIsSplit(TEST_READS.get(3),0L));
    }

    @Test(groups = "spark")
    void testDiscordantPairDetector() {
        final FindSVBreakpointsSpark.DiscordantPairDetector detector =
                new FindSVBreakpointsSpark.DiscordantPairDetector();
        Assert.assertNull(detector.getLocusIfDiscordant(TEST_READS.get(0),0L));
        Assert.assertNull(detector.getLocusIfDiscordant(TEST_READS.get(1),0L));
        Assert.assertNotNull(detector.getLocusIfDiscordant(TEST_READS.get(2),0L));
        Assert.assertNotNull(detector.getLocusIfDiscordant(TEST_READS.get(3),0L));
    }

    @Test(groups = "spark")
    void testBreakpointClusterer() {
        final FindSVBreakpointsSpark.BreakpointClusterer clusterer = new FindSVBreakpointsSpark.BreakpointClusterer();
        final GATKRead read = TEST_READS.get(0);
        // starting with 1 to do this loop one time fewer than necessary to spill reads
        for ( int nTimes = 1; nTimes != FindSVBreakpointsSpark.BreakpointClusterer.MIN_EVIDENCE; ++nTimes ) {
            Assert.assertFalse(clusterer.apply(read).hasNext());
        }
        final Iterator<GATKRead> readIterator = clusterer.apply(read);
        for ( int nTimes = 0; nTimes != FindSVBreakpointsSpark.BreakpointClusterer.MIN_EVIDENCE; ++nTimes ) {
            Assert.assertTrue(readIterator.hasNext());
            Assert.assertEquals(read, readIterator.next());
        }
        Assert.assertFalse(readIterator.hasNext());
    }

    @Test(groups = "spark")
    void testClusterKmerizer() {
        final FindSVBreakpointsSpark.ClusterKmerizer clusterKmerizer =
                new FindSVBreakpointsSpark.ClusterKmerizer(0, Collections.emptySet());
        final Iterator<Tuple2<SVKmer, Long>> iterator = clusterKmerizer.apply(TEST_READS.get(0));
        final HashSet<Tuple2<SVKmer, Long>> tupleSet = new HashSet<>();
        while ( iterator.hasNext() ) {
            tupleSet.add(iterator.next());
        }
        Assert.assertEquals(tupleSet.size(), 1);
        final Tuple2<SVKmer, Long> expectedValue = new Tuple2<>(new SVKmer(SVConstants.KMER_SIZE), 1L);
        Assert.assertTrue(tupleSet.contains(expectedValue));
    }

    @Test(groups = "spark")
    void testReadsForBreakpointFinder() {
        final Long breakpointId = 1L;
        final Map<SVKmer, List<Long>> kmerBreakpointMap = new HashMap<>();
        kmerBreakpointMap.put(new SVKmer(SVConstants.KMER_SIZE), Collections.singletonList(breakpointId));
        final FindSVBreakpointsSpark.ReadsForBreakpointFinder finder =
                new FindSVBreakpointsSpark.ReadsForBreakpointFinder(kmerBreakpointMap);
        final GATKRead read = TEST_READS.get(0);
        final Iterator<Tuple2<Long, GATKRead>> iterator = finder.apply(read);
        Assert.assertTrue(iterator.hasNext());
        final Tuple2<Long, GATKRead> expectedValue = new Tuple2<>(breakpointId, read);
        Assert.assertEquals(iterator.next(), expectedValue);
        Assert.assertFalse(iterator.hasNext());
    }
}
