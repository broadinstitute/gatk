package org.broadinstitute.hellbender.tools.exome.orientationbiasvariantfilter;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static org.mockito.Matchers.any;
import static org.mockito.Mockito.*;

public class OxoQScorerUnitTest extends BaseTest {

    public static final String largeFileAcnvFileTestDir = largeFileTestDir + "/ACNV_pipeline_test_files/";
    public static final String testBam = largeFileAcnvFileTestDir + "HCC1143-t1-chr20-downsampled.deduplicated.bam";
    public static final String testRef = largeFileAcnvFileTestDir + "human_g1k_v37.chr-20.truncated.fasta";

    /** This is testing a private method, so there is not much error testing, since calling functions will make sure
     * inputs are valid.
     *
     * Three keys should be generated, since we assume that the ref bases are padded by one base on either side and
     *  the contexts are only one base long
     */
    @Test
    public void testBasicOxoQBinKeyGeneration() throws IOException {
        final String readBases = "CGATAGGT";
        final String refBases = "ACGCTACCAA";
        final byte[] fakeReadQuals = new byte[readBases.length()];
        Arrays.fill(fakeReadQuals, (byte) 50);

        SimpleInterval interval = new SimpleInterval("22", 500, 500+refBases.length()-1);

        ReferenceMultiSource mockSource = mock(ReferenceMultiSource.class, withSettings().serializable());
        when(mockSource.getReferenceBases(any(PipelineOptions.class), any(SimpleInterval.class))).thenReturn(new ReferenceBases(refBases.getBytes(), interval));

        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

        // fake, but accurate cigar
        final String fakeCigarString = CigarUtils.calculateCigar(Arrays.copyOfRange(refBases.getBytes(), 1, refBases.length()-1), readBases.getBytes()).toString();
        GATKRead read = ArtificialReadUtils.createArtificialRead(header, "fake_read", interval.getContig(), interval.getStart() + 1, readBases.getBytes(), fakeReadQuals, fakeCigarString);

        final List<OxoQBinKey> result = OxoQScorer.createOxoQBinKeys(new Tuple2<>(read, mockSource.getReferenceBases((PipelineOptions) null, interval)), OxoQScorer.createStringOxoQBinKeyMap());
        Assert.assertNotNull(result);
        Assert.assertEquals(result.size(), read.getLength() - 2 );
        Assert.assertEquals(result.size(), readBases.length() - 2);
        for (int i = 0; i < read.getLength()-2; i ++) {
            Assert.assertEquals(result.get(i), new OxoQBinKey(new ArtifactMode(refBases.getBytes()[i+2], readBases.getBytes()[i+1]),
                    new Tuple2<>(readBases.getBytes()[i], readBases.getBytes()[i+2]), true));
        }
    }

    /** Just make sure we can create the OxoQKeys from a bam file.  Just checks for crashes and little else. */
    @Test(timeOut = 60000)
    public void testBasic() {
        final ReferenceMultiSource initialRef = new ReferenceMultiSource((PipelineOptions) null, testRef, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final ReferenceMultiSource reference = new ReferenceMultiSource((PipelineOptions) null, testRef, new OxoQScorer.OxoQBinReferenceWindowFunction(initialRef.getReferenceSequenceDictionary(null)));

        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        ReadsSparkSource src = new ReadsSparkSource(ctx);
        final ReadFilter filter = makeGenomeReadFilter();
        JavaRDD<GATKRead> rawReads = src.getParallelReads(testBam, testRef);
        JavaRDD<GATKRead> reads = rawReads.filter(read -> filter.test(read));
        final double score = OxoQScorer.scoreReads(reads, reference, ctx);
        Assert.assertTrue(score >= 0);
    }

    @Test(timeOut = 60000)
    public void testScoreInBasicScenario() throws IOException{
        final String read1Bases = "CGATAGGT";
        final String read2Bases = "CGCTAGGT";
        final String refBases =  "ACGCTACCAA";

        final byte[] fakeReadQuals = new byte[read1Bases.length()];
        Arrays.fill(fakeReadQuals, (byte) 50);

        SimpleInterval interval = new SimpleInterval("22", 500, 500+refBases.length()-1);

        ReferenceMultiSource mockSource = mock(ReferenceMultiSource.class, withSettings().serializable());
        when(mockSource.getReferenceBases(any(PipelineOptions.class), any(SimpleInterval.class))).thenReturn(new ReferenceBases(refBases.getBytes(), interval));

        SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

        // fake, but accurate cigar
        final String fakeCigarString = CigarUtils.calculateCigar(Arrays.copyOfRange(refBases.getBytes(), 1, refBases.length()-1), read1Bases.getBytes()).toString();
        GATKRead read1 = ArtificialReadUtils.createArtificialRead(header, "fake_read1", interval.getContig(), interval.getStart() + 1, read1Bases.getBytes(), fakeReadQuals, fakeCigarString);
        GATKRead read2 = ArtificialReadUtils.createArtificialRead(header, "fake_read2", interval.getContig(), interval.getStart() + 1, read2Bases.getBytes(), fakeReadQuals, fakeCigarString);
        read2.setIsSecondOfPair();

        final List<GATKRead> readsAsList = new LinkedList<>();
        readsAsList.add(read1);
        readsAsList.add(read2);

        JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
        JavaRDD<GATKRead> reads = ctx.parallelize(readsAsList);
        final double score = OxoQScorer.scoreReads(reads, mockSource, ctx);

        //TODO: Test a raw score...


    }

    private ReadFilter makeGenomeReadFilter() {
        return ReadFilterLibrary.VALID_ALIGNMENT_START
                .and(ReadFilterLibrary.VALID_ALIGNMENT_END)
                .and(ReadFilterLibrary.HAS_READ_GROUP)
                .and(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS)
                .and(ReadFilterLibrary.READLENGTH_EQUALS_CIGARLENGTH)
                .and(ReadFilterLibrary.SEQ_IS_STORED)
                .and(ReadFilterLibrary.CIGAR_CONTAINS_NO_N_OPERATOR)
                .and(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT)
                .and(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK)
                .and(ReadFilterLibrary.MAPPED);
    }
}
