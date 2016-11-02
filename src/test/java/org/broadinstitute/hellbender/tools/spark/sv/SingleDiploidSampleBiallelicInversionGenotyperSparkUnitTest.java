package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.SingleDiploidSampleBiallelicInversionGenotyperSpark.InversionJunction;

/**
 * Created by shuang on 9/2/16.
 */
public final class SingleDiploidSampleBiallelicInversionGenotyperSparkUnitTest extends BaseTest{

    private static JavaSparkContext ctx;

    private ReferenceMultiSource referenceMultiSource;

    private List<VariantContext> vcs;
    private List<InversionJunction> junctions;
    private final Map<Long, List<LocalAssemblyContig>> assemblyID2assembleContents = new HashMap<>(); // empty because currently this is not used
    private static final SingleDiploidSampleBiallelicInversionGenotyperSpark tool = new SingleDiploidSampleBiallelicInversionGenotyperSpark();

    @BeforeClass
    private void setupSparkAndTestFile(){
        SparkContextFactory.enableTestSparkContext();
        ctx = SparkContextFactory.getTestSparkContext(Collections.emptyMap());

        referenceMultiSource = new ReferenceMultiSource((PipelineOptions)null, new File(b37_reference_20_21).getAbsolutePath(), ReferenceWindowFunctions.IDENTITY_FUNCTION);

        final String inputVCF = new File("src/test/resources/org/broadinstitute/hellbender/tools/spark/sv/SingleDiploidSampleBiallelicSVGenotyperSpark").getAbsolutePath() + "/inversions.vcf";
        try(final VCFFileReader reader = new VCFFileReader(new File(inputVCF), false)){
            vcs = StreamSupport.stream(reader.spliterator(), false).collect(Collectors.toList());
            junctions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());
        }
    }

    @AfterClass
    private static void closeSpark(){
        SparkContextFactory.stopSparkContext(ctx);
    }

    @Test
    public void testReadSuitableForGenotypingJunction(){

        final InversionJunction junction = junctions.get(0);

        final String fakeSeq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        final String fakeQual = "#######################################################################################################################################################";
        final GATKRead fastqRead = SVFastqUtils.convertToRead(new FastqRecord("11259L mapping=unmapped", fakeSeq, "+", fakeQual));

        Assert.assertTrue(tool.readSuitableForGenotypingJunction(junction, fastqRead));

        final GATKRead unmappedRead = ArtificialReadUtils.createRandomRead(151); unmappedRead.setName("unmappedRead"); unmappedRead.setIsUnmapped(); unmappedRead.setIsSecondOfPair();
        Assert.assertFalse( tool.readResideInJunctionWindow(junction, unmappedRead) );

        final GATKRead wrongChr = ArtificialReadUtils.createRandomRead(151);wrongChr.setName("wrongChr"); wrongChr.setPosition("21", 27374100);wrongChr.setIsFirstOfPair();
        Assert.assertFalse( tool.readResideInJunctionWindow(junction, wrongChr) );

        final GATKRead tooFarLeft = ArtificialReadUtils.createRandomRead(151);tooFarLeft.setName("tooFarLeft");tooFarLeft.setPosition("20", 26210700); tooFarLeft.setIsFirstOfPair(); // too far left
        Assert.assertFalse( tool.readResideInJunctionWindow(junction, tooFarLeft) );

        final GATKRead tooFarRight = ArtificialReadUtils.createRandomRead(151);tooFarRight.setName("tooFarRight"); tooFarRight.setPosition("20", 26210908);tooFarRight.setIsSecondOfPair(); // too far right
        Assert.assertFalse( tool.readResideInJunctionWindow(junction, tooFarRight) );

        final GATKRead correctRead = ArtificialReadUtils.createRandomRead(151); correctRead.setName("correctRead"); correctRead.setPosition("20", 26210857); correctRead.setIsFirstOfPair();
        Assert.assertTrue( tool.readResideInJunctionWindow(junction, correctRead) );

    }

    @Test
    public void testConvertToSVJunction(){

        Assert.assertEquals(junctions.size(), 5);

        final InversionJunction junction =  junctions.get(0);

        // test only one element, not all
        Assert.assertEquals(junction.getOriginalVC(), vcs.get(0));
        Assert.assertEquals(junction.getOriginalVC(), vcs.get(0));
        Assert.assertEquals(junction.type, GATKSVVCFHeaderLines.SVTYPES.INV);
        Assert.assertEquals(junction.fiveEndBPLoc, new SimpleInterval("20", 26210907, 26210907));
        Assert.assertEquals(junction.threeEndBPLoc, new SimpleInterval("20", 26211109, 26211109));
        Assert.assertEquals(junction.assemblyIDs, Arrays.asList(5718L, 7161L));
        Assert.assertEquals(junction.contigIDs, Arrays.asList(new Tuple2<>(5718L, "contig-7"), new Tuple2<>(7161L, "contig-3")));
        Assert.assertEquals(junction.numHQMappings, 2);
        Assert.assertEquals(junction.maxAlignLength, 81);
        Assert.assertEquals(junction.svLength, 202);
        Assert.assertNull(junction.maybeNullHomology);
        Assert.assertEquals(junction.maybeNullInsertedSeq, "CACACACACACGTAGCCTCATAATACCATATATATATA".getBytes());
        Assert.assertEquals(junction.invType, BreakpointAllele.InversionType.INV_5_TO_3);
    }

    @Test
    public void testUpdateReads(){

        final InversionJunction junction = junctions.get(0);
        // set up test data
        final List<GATKRead> reads = new ArrayList<>();
        final ReadLikelihoods<SVDummyAllele> input = readPostprocessingTestDataBuilder(junction, reads);

        // because normalization, filtering and marginalization are all tested, we simply test the end output
        final ReadLikelihoods<SVDummyAllele> output = tool.updateReads(junction, input);
        Assert.assertEquals(output.samples().size(), 1);
        Assert.assertEquals(output.samples().get(0), SingleDiploidSampleBiallelicSVGenotyperSpark.testSampleName);

        final LikelihoodMatrix<SVDummyAllele> mat = output.sampleMatrix(0);
        Assert.assertEquals(mat.alleles(), junction.getGenotypedVC().getAlleles()); // test alleles

        Assert.assertEquals(mat.reads().size(), 4); // no reads filtered out because test data is designed to be well modelled
        final List<Double> firstRow = Arrays.asList(mat.get(0, 0), mat.get(0, 1), mat.get(0, 2), mat.get(0, 3));
        Assert.assertTrue(firstRow.indexOf(Collections.max(firstRow))<2);
        final List<Double> secondRow = Arrays.asList(mat.get(1, 0), mat.get(1, 1), mat.get(1, 2), mat.get(1, 3));
        Assert.assertTrue(secondRow.indexOf(Collections.max(secondRow))>1);
        Assert.assertTrue(mat.get(0, 2) < -1.0);
        Assert.assertTrue(mat.get(0, 3) < -1.0);
        Assert.assertTrue(mat.get(1, 0) < -1.0);
        Assert.assertTrue(mat.get(1, 1) < -1.0);
    }

    // set up test data
    private static ReadLikelihoods<SVDummyAllele> readPostprocessingTestDataBuilder(final InversionJunction junction, final List<GATKRead> reads) {

        reads.addAll(IntStream.range(1, 5).mapToObj(i -> {
            final GATKRead read = ArtificialReadUtils.createRandomRead(151); read.setName("read_"+i); return read;
        }).collect(Collectors.toList()));

        final Map<String, List<GATKRead>> sample2Reads = new HashMap<>();
        sample2Reads.put(SingleDiploidSampleBiallelicSVGenotyperSpark.testSampleName, reads);

        final ReadLikelihoods<SVDummyAllele> rll = new ReadLikelihoods<>(new IndexedSampleList(SingleDiploidSampleBiallelicSVGenotyperSpark.testSampleName),
                new IndexedAlleleList<>(junction.getAlleles()), sample2Reads);

        final int alleleCnt = junction.getAlleles().size();

        // copied from ReadLikelihoodsUnitTest.java with modification of Allele type
        Utils.resetRandomGenerator();
        final Random rnd = Utils.getRandomGenerator();
        final double[][] likelihoods = new double[alleleCnt][reads.size()];
        final LikelihoodMatrix<SVDummyAllele> sampleLikelihoods = rll.sampleMatrix(0);
        for (int a = 0; a < alleleCnt; a++) {
            for (int r = 0; r < likelihoods[0].length; r++)
                sampleLikelihoods.set(a,r,likelihoods[a][r] = -1.0-Math.abs(rnd.nextGaussian()));
        }
        sampleLikelihoods.set(0, 0, likelihoods[0][0] = -0.1 + rnd.nextGaussian()*0.01);
        sampleLikelihoods.set(1, 1, likelihoods[1][1] = -0.1 + rnd.nextGaussian()*0.01);
        sampleLikelihoods.set(2, 2, likelihoods[2][2] = -0.1 + rnd.nextGaussian()*0.01);
        sampleLikelihoods.set(3, 3, likelihoods[3][3] = -0.1 + rnd.nextGaussian()*0.01);

        return rll;
    }
}
