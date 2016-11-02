package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.lang3.StringUtils;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.InversionType;
import static org.broadinstitute.hellbender.tools.spark.sv.SingleDiploidSampleBiallelicInversionGenotyperSpark.InversionJunction;
/**
 * Created by shuang on 9/16/16.
 */
public final class SVJunctionUnitTest extends BaseTest{

    private static JavaSparkContext ctx;

    private ReferenceMultiSource referenceMultiSource;

    private List<VariantContext> vcs;
    private List<InversionJunction> inversions;

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
            inversions = vcs.stream().map(vc -> tool.convertToSVJunction(vc, ctx.broadcast(assemblyID2assembleContents), ctx.broadcast(referenceMultiSource))).collect(Collectors.toList());
        }

    }

    @Test
    public void testGetAssociatedContigs(){
        final AlignmentRegion expected = new AlignmentRegion("1", "contig-0", TextCigarCodec.decode("151M"), true, new SimpleInterval("20", 1000000, 1000000), 60, 1, 151, 0);
        final LocalAssemblyContig onezero = new LocalAssemblyContig(1L, "contig-0", "AAA", new ArrayList<>(Collections.singletonList(expected)));
        final LocalAssemblyContig oneone = new LocalAssemblyContig(1L, "contig-1", "CCC");
        final LocalAssemblyContig onetwo = new LocalAssemblyContig(1L, "contig-2", "GGG");
        final LocalAssemblyContig onethree = new LocalAssemblyContig(1L, "contig-3", "TTT");

        final LocalAssemblyContig twozero = new LocalAssemblyContig(2L, "contig-0", "TTT");
        final LocalAssemblyContig twoone = new LocalAssemblyContig(2L, "contig-1", "GGG");
        final LocalAssemblyContig twotwo = new LocalAssemblyContig(2L, "contig-2", "CCC");
        final LocalAssemblyContig twothree = new LocalAssemblyContig(2L, "contig-3", "AAA");

        // first construct input map
        final Map<Long, List<LocalAssemblyContig>> inputMap = new HashMap<>();
        inputMap.put(1L, Arrays.asList(onezero, oneone, onetwo, onethree));
        inputMap.put(2L, Arrays.asList(twozero, twoone, twotwo, twothree)); // opposite order of 1L

        // second construct input list
        final List<Tuple2<Long, String>> inputList = Arrays.asList(new Tuple2<>(1L, "contig-0"), new Tuple2<>(2L, "contig-1"), new Tuple2<>(2L, "contig-2"), new Tuple2<>(2L, "contig-0"));

        final List<LocalAssemblyContig> tobeTested = SVJunction.getAssociatedContigs(inputMap, inputList);

        Assert.assertEquals(tobeTested.size(), 4);

        Assert.assertTrue(tobeTested.containsAll(Arrays.asList(onezero, twozero, twoone, twotwo)));
    }

    //////////////////////// sub classes tests

    @Test
    public void testInversionEquals(){
        final InversionJunction invOne = inversions.get(0);
        final InversionJunction invTwo = inversions.get(1);
        Assert.assertNotEquals(invOne, invTwo);
    }

    @Test
    public void testInversionHashcode(){
        final InversionJunction invOne = inversions.get(0);
        final InversionJunction invTwo = inversions.get(1);
        Assert.assertNotEquals(invOne.hashCode(), invTwo.hashCode());
    }

    @Test
    public void testInversionSetEnd(){
        Assert.assertEquals(inversions.stream().map(j -> j.whichEnd).toArray(), new byte[]{1, 1, 1, 1, 2});
    }

    @Test
    public void testInversionConstructReferenceWindows(){
        final int readL = SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength;

        inversions.forEach( inv -> {
            int flank = inv.maybeNullHomology == null ? 0 : inv.maybeNullHomology.length;

            Tuple2<SimpleInterval, SimpleInterval> refWindows = inv.getReferenceWindows();
            Assert.assertEquals(refWindows._1().getContig(), refWindows._2().getContig());
            Assert.assertEquals(inv.getLeftAlignedBreakpointLocations().get(0).getStart() - refWindows._1().getStart(), readL-1);
            Assert.assertEquals(refWindows._1().getEnd() - inv.getLeftAlignedBreakpointLocations().get(0).getStart(), readL + flank * (inv.invType== InversionType.INV_NONE ? 0 : 1) - 2);
            Assert.assertEquals(inv.getLeftAlignedBreakpointLocations().get(1).getStart() - refWindows._2().getStart(), readL + flank * (inv.invType== InversionType.INV_NONE ? 0 : 1) - 1);
            Assert.assertEquals(refWindows._2().getEnd() - inv.getLeftAlignedBreakpointLocations().get(1).getStart(), readL-2);
        });
    }

    @Test
    public void testInversionConstructRefAlleles(){

        try (final ContigAligner contigAligner = new ContigAligner(b37_reference_20_21)) {
            inversions.forEach( inv -> {
                final List<SVDummyAllele> alleles = inv.getAlleles();

                testRefAllelesByAlignment(alleles.get(0), true, inv, contigAligner);
                testRefAllelesByAlignment(alleles.get(1), false, inv, contigAligner);
                Assert.assertEquals(alleles.get(0).length(), alleles.get(1).length());
            });
        }catch (final IOException e) {
            throw new GATKException("Cannot run BWA-MEM", e);
        }
    }

    private void testRefAllelesByAlignment(final SVDummyAllele allele, final boolean is5Side, final InversionJunction inv, final ContigAligner contigAligner){
        try{
            final int offsetTobeTurnedOffWhenJBWAIsFixed = 1;
            final byte[] bases = allele.getBases();
            final byte[] qual = new byte[bases.length];
            final ShortRead alleleRead = new ShortRead(inv.getOriginalVC().getID(), bases, qual);
            final AlnRgn[] alnRgns = contigAligner.bwaMem.align(alleleRead);
            Assert.assertTrue(alnRgns.length!=0);
            Assert.assertEquals(alnRgns[0].getPos()+offsetTobeTurnedOffWhenJBWAIsFixed, is5Side ? inv.getReferenceWindows()._1().getStart() : inv.getReferenceWindows()._2().getStart() );
            Assert.assertEquals(alnRgns[0].getCigar(), String.valueOf(bases.length)+"M");
        } catch (final IOException e) {
            throw new GATKException("Cannot run BWA-MEM", e);
        }
    }

    @Test
    public void testInversionConstructAltAlleles(){
        inversions.forEach( inv -> {
            final List<SVDummyAllele> alleles = inv.getAlleles();
            final int ins = inv.maybeNullInsertedSeq==null?  0 : inv.maybeNullInsertedSeq.length;
            testAltAlleles(inv, alleles.get(2), true);
            testAltAlleles(inv, alleles.get(3), false);
            Assert.assertEquals(alleles.get(2).length()-ins, alleles.get(0).length());
        });
    }

    private void testAltAlleles(final InversionJunction inv, final SVDummyAllele allele, final boolean is5Side){

        final int hom = inv.maybeNullHomology==null ? 0 : inv.maybeNullHomology.length;
        final int ins = inv.maybeNullInsertedSeq==null ? 0 : inv.maybeNullInsertedSeq.length;
        Assert.assertEquals(allele.length(), SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength*2+ins+hom-2);

        if(is5Side){
            // test first "half" of the alt allele is the same as the first "half" of the left ref allele
            byte[] s = Arrays.copyOfRange(inv.getAlleles().get(0).getBases(), 0, SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1);
            final int distL = StringUtils.getLevenshteinDistance(new String(Arrays.copyOfRange(allele.getBases(), 0, SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1)),
                                                                 new String(s));
            Assert.assertEquals(distL, 0);

            // test second "half" of the alt allele is the same as RC of the first "half" of the right ref allele
                   s = Arrays.copyOfRange(inv.getAlleles().get(1).getBases(), 0, SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1);
            SequenceUtil.reverseComplement(s);
            final int distR = StringUtils.getLevenshteinDistance(new String(Arrays.copyOfRange(allele.getBases(), SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1 + hom + ins, allele.length())), // skip insertion and homology
                                                                 new String(s));
            Assert.assertEquals(distR, 0);
        } else {
            // test first "half" of the alt allele is the same as the second "half" of the left ref allele
            SVDummyAllele a = inv.getAlleles().get(0);
            byte[] s = Arrays.copyOfRange(a.getBases(), SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1+hom, a.length());                                                                        // skip homology in ref
            SequenceUtil.reverseComplement(s);
            final int distL = StringUtils.getLevenshteinDistance(new String(Arrays.copyOfRange(allele.getBases(), 0, SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1)),
                                                                 new String(s));
            Assert.assertEquals(distL, 0);

            // test second "half" of the alt allele is the same as the second "half" of the right ref allele
            a = inv.getAlleles().get(1);
            s = Arrays.copyOfRange(a.getBases(), SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1+hom, a.length());
            final int distR = StringUtils.getLevenshteinDistance(new String(Arrays.copyOfRange(allele.getBases(), SingleDiploidSampleBiallelicInversionGenotyperSpark.readLength-1 + hom + ins, allele.length())), // skip insertion and homology
                                                                 new String(s));
            Assert.assertEquals(distR, 0);
        }
    }
}
