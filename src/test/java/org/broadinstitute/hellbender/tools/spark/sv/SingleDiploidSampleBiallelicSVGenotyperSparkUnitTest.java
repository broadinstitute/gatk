package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.io.FileUtils;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Created by shuang on 9/2/16.
 */
public final class SingleDiploidSampleBiallelicSVGenotyperSparkUnitTest extends BaseTest{

    private JavaSparkContext ctx;

    private ReferenceMultiSource referenceMultiSource;

    private List<VariantContext> vcs;

    private final Map<Long, List<LocalAssemblyContig>> assemblyID2assembleContents = new HashMap<>(); // empty because currently this is not used
    private static final SingleDiploidSampleBiallelicInversionGenotyperSpark tool = new SingleDiploidSampleBiallelicInversionGenotyperSpark();

    @BeforeClass
    private void setupSparkAndTestFile(){
        SparkContextFactory.enableTestSparkContext();
        ctx = SparkContextFactory.getTestSparkContext(Collections.emptyMap());

        referenceMultiSource = new ReferenceMultiSource((PipelineOptions)null, new File(b37_reference_20_21).getAbsolutePath(), ReferenceWindowFunctions.IDENTITY_FUNCTION);

        final String inputVCF = new File(getToolTestDataDir()).getAbsolutePath() + "/inversions.vcf";
        try(final VCFFileReader reader = new VCFFileReader(new File(inputVCF), false)){
            vcs = StreamSupport.stream(reader.spliterator(), false).collect(Collectors.toList());
        }
    }

    @AfterClass
    private void closeSpark(){
        SparkContextFactory.stopSparkContext(ctx);
    }

    @Test
    public void testGatherAsmIDs(){
        final Set<Long> asmIDs = SingleDiploidSampleBiallelicSVGenotyperSpark.gatherAsmIDs(vcs);
        Assert.assertEquals(asmIDs.size(), 7);
        Assert.assertTrue(asmIDs.containsAll(Arrays.asList(5718L, 7161L, 11259L, 12583L, 12697L, 12944L, 13670L)));
    }

    @Test
    public void testMapAssemblyID2ItsAlignments() throws IOException{

        // assembly
        final File testFASTAFile = new File("src/test/resources/org/broadinstitute/hellbender/tools/spark/sv/RunSGAViaProcessBuilderOnSpark/assembly9.pp.ec.filter.pass.merged.rmdup-contigs.fa");
        final ContigsCollection collection = new ContigsCollection(1, FileUtils.readLines(testFASTAFile, "UTF-8"));
        final File tempFile = createTempFile("whatever", ""); tempFile.deleteOnExit();
        FileUtils.writeStringToFile(tempFile, "9\t"+collection.toPackedFasta());

        final JavaPairRDD<Long, ContigsCollection> asmId2Contigs = ContigsCollection.loadContigsCollectionKeyedByAssemblyId(ctx, "file://"+tempFile.getAbsolutePath())
                .mapToPair(pair -> new Tuple2<>(Long.parseLong(pair._1().replace(SVConstants.FASTQ_OUT_PREFIX, "")), pair._2()));

        // alignments for this file
        final AlignmentRegion contig4Ar1 = ContigsCollection.parseAlignedAssembledContigLine("9\tcontig-4\tX\t14731120\t14731319\t+\t200M149S\t60\t1\t200\t0");
        final AlignmentRegion contig4Ar2 = ContigsCollection.parseAlignedAssembledContigLine("9\tcontig-4\tX\t14729392\t14729541\t-\t199H150M\t60\t200\t349\t0");
        final AlignmentRegion contig1Ar1 = ContigsCollection.parseAlignedAssembledContigLine("9\tcontig-1\tX\t14731320\t14731464\t-\t145M143S\t60\t1\t145\t0");
        final AlignmentRegion contig1Ar2 = ContigsCollection.parseAlignedAssembledContigLine("9\tcontig-1\tX\t14729539\t14729682\t+\t144H144M\t60\t145\t288\t0");

        List<Tuple2<Long, Iterable<AlignmentRegion>>> alignments = Arrays.asList(new Tuple2<>(9L, Arrays.asList(contig4Ar1, contig4Ar2, contig1Ar1, contig1Ar2)));
        JavaPairRDD<Long, Iterable<AlignmentRegion>> asmId2Alignments = ctx.parallelizePairs(alignments);

        final Map<Long, List<LocalAssemblyContig>> result = SingleDiploidSampleBiallelicSVGenotyperSpark.mapAssemblyID2ItsAlignments(asmId2Contigs, asmId2Alignments, Collections.singleton(9L));

        Assert.assertEquals(result.size(), 1);
        Assert.assertEquals(result.keySet().iterator().next(), Long.valueOf(9));
        final List<LocalAssemblyContig> contigs = result.get(9L);
        Assert.assertEquals(contigs.size(), 7);
        Collections.sort(contigs, (LocalAssemblyContig contig1, LocalAssemblyContig contig2) -> contig1.contigID.compareTo(contig2.contigID));
        Assert.assertNull(contigs.get(0).alignmentRecords);
        Assert.assertNull(contigs.get(2).alignmentRecords);
        Assert.assertNull(contigs.get(3).alignmentRecords);
        Assert.assertNull(contigs.get(5).alignmentRecords);
        Assert.assertNull(contigs.get(6).alignmentRecords);

        Assert.assertTrue(contigs.get(1).alignmentRecords.contains(contig1Ar1));
        Assert.assertTrue(contigs.get(1).alignmentRecords.contains(contig1Ar2));
        Assert.assertTrue(contigs.get(4).alignmentRecords.contains(contig4Ar1));
        Assert.assertTrue(contigs.get(4).alignmentRecords.contains(contig4Ar2));
    }

    @Test
    public void testFilterJunctionsBasedOnIntervalKillList(){

        final List<SVJunction> juntions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());
        final Predicate<SVJunction> filter = new SingleDiploidSampleBiallelicSVGenotyperSpark.JunctionFilterBasedOnInterval();

        SingleDiploidSampleBiallelicSVGenotyperSpark.filterJunctions(juntions, filter);
        Assert.assertEquals(juntions.size(), 5);
        final Predicate<SVJunction> anotherFilter = new SingleDiploidSampleBiallelicSVGenotyperSpark.JunctionFilterBasedOnInterval(Arrays.asList(new SimpleInterval("20", 26211000, 26211000)));
        SingleDiploidSampleBiallelicSVGenotyperSpark.filterJunctions(juntions, anotherFilter);
        Assert.assertEquals(juntions.size(), 3);
    }

    @Test
    public void testGetAssemblyToJunctionsMapping(){

        final List<SVJunction> junctions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());
        final Map<Long, List<SVJunction>> assemblyToJunctionsMapping = SingleDiploidSampleBiallelicSVGenotyperSpark.getAssemblyToJunctionsMapping(junctions);

        Assert.assertEquals(assemblyToJunctionsMapping.size(), 7);
        Assert.assertTrue(assemblyToJunctionsMapping.keySet().containsAll(Arrays.asList(5718L, 7161L, 11259L, 12583L, 12697L, 12944L, 13670L)));

        Assert.assertEquals(assemblyToJunctionsMapping.get(5718L).size(), 2);
        Assert.assertEquals(assemblyToJunctionsMapping.get(5718L).get(0), junctions.get(0));
        Assert.assertEquals(assemblyToJunctionsMapping.get(5718L).get(1), junctions.get(1));

        Assert.assertEquals(assemblyToJunctionsMapping.get(7161L).size(), 1);
        Assert.assertEquals(assemblyToJunctionsMapping.get(7161L).get(0), junctions.get(0));

        Assert.assertEquals(assemblyToJunctionsMapping.get(11259L).size(), 1);
        Assert.assertEquals(assemblyToJunctionsMapping.get(11259L).get(0), junctions.get(1));

        Assert.assertEquals(assemblyToJunctionsMapping.get(12583L).size(), 1);
        Assert.assertEquals(assemblyToJunctionsMapping.get(12583L).get(0), junctions.get(1));

        Assert.assertEquals(assemblyToJunctionsMapping.get(13670L).size(), 1);
        Assert.assertEquals(assemblyToJunctionsMapping.get(13670L).get(0), junctions.get(1));

        Assert.assertEquals(assemblyToJunctionsMapping.get(12697L).size(), 1);
        Assert.assertEquals(assemblyToJunctionsMapping.get(12697L).get(0), junctions.get(2));

        Assert.assertEquals(assemblyToJunctionsMapping.get(12944L).size(), 2);
        Assert.assertEquals(assemblyToJunctionsMapping.get(12944L).get(1), junctions.get(4));
    }

    @Test
    public void testGetFASTQReads(){

        final String fakeSeq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        final String fakeQual = "#######################################################################################################################################################";

        final FastqRecord fastqFor5718L = new FastqRecord("5718L mapping=20:2621857;151M", fakeSeq, "+", fakeQual);
        final FastqRecord fastqFor7161L = new FastqRecord("7161L mapping=unmapped", fakeSeq, "+", fakeQual); // mapped or not, FASTQ reads should not be filtered out
        final FastqRecord fastqFor11259L = new FastqRecord("11259L mapping=20:2621807;151M", fakeSeq, "+", fakeQual);
        final FastqRecord fastqFor12583L = new FastqRecord("12583L mapping=unmapped", fakeSeq, "+", fakeQual);
        final FastqRecord fastqFor13670L = new FastqRecord("13670L mapping=unmapped", fakeSeq, "+", fakeQual);
        final FastqRecord fastqFor12697L = new FastqRecord("12697L mapping=20:48450800;151M", fakeSeq, "+", fakeQual);
        final FastqRecord fastqFor12944L = new FastqRecord("12944L mapping=21:27374100;151M", fakeSeq, "+", fakeQual);

        final List<Tuple2<Long, String>> asmToFASTQ = Arrays.asList(new Tuple2<>(5718L, fastqFor5718L.toString()), new Tuple2<>(7161L, fastqFor7161L.toString()),
                new Tuple2<>(11259L, fastqFor11259L.toString()), new Tuple2<>(12583L, fastqFor12583L.toString()),
                new Tuple2<>(13670L, fastqFor13670L.toString()), new Tuple2<>(12697L, fastqFor12697L.toString()), new Tuple2<>(12944L, fastqFor12944L.toString()));

        final JavaPairRDD<Long, String> fastq = ctx.parallelizePairs(asmToFASTQ);
        final List<SVJunction> svJunctions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());
        final JavaPairRDD<Long, List<SVJunction>> asmid2SVJunctions =
                ctx.parallelizePairs(SingleDiploidSampleBiallelicSVGenotyperSpark.getAssemblyToJunctionsMapping(svJunctions).entrySet().stream().map(e -> new Tuple2<>(e.getKey(), e.getValue())).collect(Collectors.toList()));

        Map<SVJunction, List<GATKRead>> fastqReads = tool.getFASTQReads(fastq, asmid2SVJunctions).collectAsMap();
        Assert.assertEquals(fastqReads.get(svJunctions.get(0)).size(), 2);
        Assert.assertTrue(fastqReads.get(svJunctions.get(0)).stream().map(GATKRead::getName).collect(Collectors.toSet()).containsAll(Arrays.asList("5718L", "7161L")));
        Assert.assertEquals(fastqReads.get(svJunctions.get(1)).size(), 4);
        Assert.assertTrue(fastqReads.get(svJunctions.get(1)).stream().map(GATKRead::getName).collect(Collectors.toSet()).containsAll(Arrays.asList("5718L", "11259L", "12583L", "13670L")));
        Assert.assertEquals(fastqReads.get(svJunctions.get(2)).size(), 1);
        Assert.assertTrue(fastqReads.get(svJunctions.get(2)).stream().map(GATKRead::getName).collect(Collectors.toSet()).containsAll(Arrays.asList("12697L")));
        Assert.assertEquals(fastqReads.get(svJunctions.get(3)).size(), 1);
        Assert.assertTrue(fastqReads.get(svJunctions.get(3)).stream().map(GATKRead::getName).collect(Collectors.toSet()).containsAll(Arrays.asList("12944L")));
        Assert.assertEquals(fastqReads.get(svJunctions.get(4)).size(), 1);
        Assert.assertTrue(fastqReads.get(svJunctions.get(4)).stream().map(GATKRead::getName).collect(Collectors.toSet()).containsAll(Arrays.asList("12944L")));
    }

    @Test
    public void testGetBamReads(){

        final GATKRead readFor5718L = ArtificialReadUtils.createRandomRead(151); readFor5718L.setName("5718L"); readFor5718L.setPosition("20", 26210857); readFor5718L.setIsFirstOfPair();
        final GATKRead readFor7161L = ArtificialReadUtils.createRandomRead(151); readFor7161L.setName("7161L"); readFor7161L.setIsUnmapped(); readFor7161L.setIsSecondOfPair();
        final GATKRead readFor11259L = ArtificialReadUtils.createRandomRead(151);readFor11259L.setName("11259L");readFor11259L.setPosition("X", 26210807); readFor11259L.setIsFirstOfPair(); // specifically mapped to wrong place
        final GATKRead readFor12583L = ArtificialReadUtils.createRandomRead(151);readFor12583L.setName("12583L"); readFor12583L.setIsUnmapped(); readFor12583L.setIsSecondOfPair();
        final GATKRead readFor13670L = ArtificialReadUtils.createRandomRead(151);readFor13670L.setName("13670L"); readFor13670L.setIsUnmapped(); readFor13670L.setIsFirstOfPair();
        final GATKRead readFor12697L = ArtificialReadUtils.createRandomRead(151);readFor12697L.setName("12697L"); readFor12697L.setPosition("20", 48450800);readFor12697L.setIsSecondOfPair();
        final GATKRead readFor12944L = ArtificialReadUtils.createRandomRead(151);readFor12944L.setName("12944L"); readFor12944L.setPosition("21", 27374100);readFor12944L.setIsFirstOfPair();

        final JavaRDD<GATKRead> reads = ctx.parallelize(new ArrayList<>(Arrays.asList(readFor5718L, readFor7161L, readFor11259L, readFor12583L, readFor13670L, readFor12697L, readFor12944L)));
        final List<SVJunction> svJunctions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());

        final Map<SVJunction, List<GATKRead>> result = tool.getBamReads(reads, svJunctions).collectAsMap();

        Assert.assertEquals(result.size(), 5);

        readFor5718L.setName("5718L/1");readFor7161L.setName("7161L/2");
        readFor11259L.setName("11259L/1");readFor12583L.setName("12583L/2");
        readFor13670L.setName("13670L/1");readFor12697L.setName("12697L/2");
        readFor12944L.setName("12944L/1");

        Assert.assertEquals(result.get(svJunctions.get(0)), Collections.singletonList(readFor5718L));
        Assert.assertEquals(result.get(svJunctions.get(1)), Collections.singletonList(readFor5718L));
        Assert.assertEquals(result.get(svJunctions.get(2)), Collections.singletonList(readFor12697L));
        Assert.assertEquals(result.get(svJunctions.get(3)), Collections.singletonList(readFor12944L));
        Assert.assertEquals(result.get(svJunctions.get(4)), Collections.singletonList(readFor12944L));
    }

    @Test
    public void testReadsMergingByName(){

        final GATKRead readOne = ArtificialReadUtils.createRandomRead(151); readOne.setName("H25T3CCXX150306:2:1105:32658:64598/1");
        final GATKRead readTwo = ArtificialReadUtils.createRandomRead(151); readTwo.setName("H25T3CCXX150306:2:1105:32658:64598/2");
        final GATKRead readThree = ArtificialReadUtils.createRandomRead(151);readThree.setName("H25T3CCXX150306:2:2124:25400:45136/1");
        final List<GATKRead> listOne = Arrays.asList(readOne, readTwo);
        final List<GATKRead> listTwo = Arrays.asList(readOne, readThree);
        final List<GATKRead> merged = SingleDiploidSampleBiallelicSVGenotyperSpark.mergeReadsByName(listOne, listTwo);
        Assert.assertEquals(merged.size(), 3);
        Assert.assertTrue(merged.containsAll(Arrays.asList(readOne, readTwo, readThree)));
    }

    @Test
    public void testInferGenotypeFromPL(){

        final List<Allele> dummyAlleles = Arrays.asList(Allele.create("A", true), Allele.create("T", false));

        GenotypeBuilder builder = new GenotypeBuilder("dummy", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(new int[] {0, 0, 0});
        Genotype gt = SingleDiploidSampleBiallelicSVGenotyperSpark.inferGenotypeFromPL(builder.make(), dummyAlleles);
        Assert.assertFalse(gt.hasPL());

        builder = new GenotypeBuilder("dummy", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(new int[] {0, 30, 500});
        gt = SingleDiploidSampleBiallelicSVGenotyperSpark.inferGenotypeFromPL(builder.make(), dummyAlleles);
        Assert.assertTrue(gt.hasPL() && gt.isHomRef());

        builder = new GenotypeBuilder("dummy", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(new int[] {30, 0, 500});
        gt = SingleDiploidSampleBiallelicSVGenotyperSpark.inferGenotypeFromPL(builder.make(), dummyAlleles);
        Assert.assertTrue(gt.hasPL() && gt.isHet());

        builder = new GenotypeBuilder("dummy", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).PL(new int[] {30, 500, 0});
        gt = SingleDiploidSampleBiallelicSVGenotyperSpark.inferGenotypeFromPL(builder.make(), dummyAlleles);
        Assert.assertTrue(gt.hasPL() && gt.isHomVar());
    }

    @Test
    public void testCollectAndSortGenotypedVC(){
        final List<SVJunction> junctions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());
        Collections.shuffle(junctions);
        final GenotypeLikelihoods justEnough = GenotypeLikelihoods.fromPLs(new int[] {0, 60, 100});
        final JavaPairRDD<SVJunction, GenotypeLikelihoods> input = ctx.parallelizePairs( junctions.stream().map(j -> new Tuple2<>(j, justEnough)).collect(Collectors.toList()) );
        final List<VariantContext> sortedGenotypedVCs = SingleDiploidSampleBiallelicSVGenotyperSpark.collectAndSortGenotypedVC(input, referenceMultiSource.getReferenceSequenceDictionary(null));
        Assert.assertEquals(sortedGenotypedVCs.size(), 5);

        Assert.assertEquals(sortedGenotypedVCs.get(0).getID(), "INV_5_TO_3_20_26210907_26211109");
        Assert.assertEquals(sortedGenotypedVCs.get(1).getID(), "INV_5_TO_3_20_26210907_26211109");
        Assert.assertEquals(sortedGenotypedVCs.get(2).getID(), "INV_5_TO_3_20_48450848_48450921");
        Assert.assertEquals(sortedGenotypedVCs.get(3).getID(), "INV_5_TO_3_21_27374151_27374706");
        Assert.assertEquals(sortedGenotypedVCs.get(4).getID(), "INV_3_TO_5_21_27374159_27374701");
    }

    @Test
    public void testOutput() throws IOException{

        final File tempGenotypedVCF = createTempFile("genotyped", "vcf"); tempGenotypedVCF.deleteOnExit();

        final List<SVJunction> junctions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());
        final int[] pls = new int[] {0, 60, 100};
        final GenotypeLikelihoods justEnough = GenotypeLikelihoods.fromPLs(pls);
        final JavaPairRDD<SVJunction, GenotypeLikelihoods> input = ctx.parallelizePairs( junctions.stream().map(j -> new Tuple2<>(j, justEnough)).collect(Collectors.toList()) );
        final List<VariantContext> sortedGenotypedVCs = SingleDiploidSampleBiallelicSVGenotyperSpark.collectAndSortGenotypedVC(input, referenceMultiSource.getReferenceSequenceDictionary(null));

        try(final OutputStream outputStream = new FileOutputStream(tempGenotypedVCF)){
            SingleDiploidSampleBiallelicSVGenotyperSpark.output(sortedGenotypedVCs, referenceMultiSource.getReferenceSequenceDictionary(null), outputStream);
        }

        try(final VCFFileReader reader = new VCFFileReader(tempGenotypedVCF, false)){
            final List<VariantContext> vcsFromFile = StreamSupport.stream(reader.spliterator(), false).collect(Collectors.toList());
            Assert.assertEquals(vcsFromFile.get(0).getID(), vcs.get(0).getID());
            Assert.assertEquals(vcsFromFile.get(0).getGenotype(0).getPL(), pls);

            Assert.assertEquals(vcsFromFile.get(1).getID(), vcs.get(1).getID());
            Assert.assertEquals(vcsFromFile.get(1).getGenotype(0).getPL(), pls);

            Assert.assertEquals(vcsFromFile.get(2).getID(), vcs.get(2).getID());
            Assert.assertEquals(vcsFromFile.get(2).getGenotype(0).getPL(), pls);

            Assert.assertEquals(vcsFromFile.get(3).getID(), vcs.get(3).getID());
            Assert.assertEquals(vcsFromFile.get(3).getGenotype(0).getPL(), pls);

            Assert.assertEquals(vcsFromFile.get(4).getID(), vcs.get(4).getID());
            Assert.assertEquals(vcsFromFile.get(4).getGenotype(0).getPL(), new int[] {0, 60, 100});
        }
    }
}
