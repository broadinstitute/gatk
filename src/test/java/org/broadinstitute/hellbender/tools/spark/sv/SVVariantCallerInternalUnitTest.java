package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceWindowFunctions;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;
import scala.Tuple3;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.BreakpointAllele.BreakpointAlleleInversion;

public class SVVariantCallerInternalUnitTest extends BaseTest {

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegions() throws Exception {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("532M87S"), true, new SimpleInterval("8", 118873207, 118873739), 60, 1, 532, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("518S29M72S"), false, new SimpleInterval("1", 175705642, 175705671), 3, 519, 547, 0);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("543S76M"), false, new SimpleInterval("1", 118875262, 118875338), 60, 544, 619, 0);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2, region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getChimericAlignmentsFromAlignmentRegions(contigSequence, alignmentRegionList, 50);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.region1, region1);
        Assert.assertEquals(chimericAlignment.region2, region3);
        Assert.assertEquals(chimericAlignment.homology, "");
        Assert.assertEquals(chimericAlignment.insertedSequence, "GAGATAGAGTC");
    }

    @Test
    public void testGetAssembledBreakpointsFromAlignmentRegionsWithOverlappingAlignmentRegion() throws Exception {
        final byte[] contigSequence = "ACTAGAGCATCTACGTGTTCCTGTGGTTTTGGAGCAAGAGTGATTTGAGTTTCAGAGATTTTTACTAATTCTTCTTCCCCTACCAGAAAAAAAGATCTTACCATTTGAGAGTGAGATGTAAACCCAGCCCTGTCTGACCTGAGTCTGTGCCCTAAGCCTATGCTAAGCCAAGCAGTGCCTGGAGCCACCACAGGTCCACACAATTCGTTAACATGATGAAGCAAGGATGGAAATTGGACAAAATAGTGTGCCTACTGAATCTAAGAATGAAAAATGATTGCACTCCTACTCTGAGTGCTTTGGAGCACTGCCCAGTTGGGCAAAGGGTCAGCGCCTGGGCAGAGGTCCCCACAACCTGGCAGGAGTGTGGTCGGCCACCCTATGGGCCTCCATCATGTGCAGTGACAGCGGGGCTGTCATGTCACCGTGTGGGAGGGCTTGCAGGTGAAGTGGTCTGGGAGGGGTCCCCCAGACAAAGCCAAGGTTCTGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCACCTGCAAGCCCTCCCACACGGTGACATGACAGCCTATAATACAGTTCCGCATGGCCACGTCATACAACCCTGTCATATTGGTGAGCAATTGCTGTGTAGCCAAAGACCCCAAAACTCAAACAGCATTTATTATTATTGCCCCCATGTCTGAGAGTCAGATGTGCATTTGCTGATCTCAGCTTGTTTGAGCTGCTGCAGGGTTGGGGCTCTGCTCCAGGCAGGCTTAGCTGTCACCACATGCACACATACATTCTGGGCCTCTGCTGCGCGCGTCACGTTCACTGAAGATCTTGGGATTGGGAGTTAGGGCGGTGGGAGGGCCCAGCAAAGTCACCTGGCGATGGCAGGGACACAGGGAGGAATGTAGAATGGGGCCGATGATGGGACCCACACGTCTGCAAAGCTGCGGTCTCCTTGAGGGGTGGAGACAGCAACAACTCACCGCACGCGGTGCTTCAGTTCACCATCTCCCTGGGACATTAGGGGGCCCCGTGTTATCTCATTTTGCTCTGGTTTGCATTAGTTTTTTATCACTTCGTAGATGAAGCCACTGACACCCAGAGAGGGAAAGTGGCCTGACCAAGGGCCACAGCAGGGGAGCGAAGGAGCCCCACAGTTCGGCAGGAACACAGCCTCTCCCTGGCTTTCAGGTTCACTGACATCTTCTCATGGCCTCTGTAACTCACCAGGCATCAGGGTGTAGTCCTTAGACCAGTGTCCCACAGCTGCCACAGAGTGGGAGCTCACCATCAGTTATAAGTCACTAGAAAGGCTTTTGGACATTATAAGCTACAATGGAAAATAAGTCATCTGTGGATTTTTGTGACAGATTCCAAAAATTTGAATATTTTGTCTACTTAGGTTTTTGGTTAATTTTATCCTCAAAACTGTTCTGCAGTGATTAAGCTGTACAAACTGCATCATGGGCGAATTGGCATATTCAGAAATGACTGATATTCTTGATTTCAGTTTTTTACTTTGTATGTAGCTCCTCAAGGAAAC".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("487M1006S"), true, new SimpleInterval("20", 23102817, 23103304), 60, 1, 487, 1);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("483S42M968S"), false, new SimpleInterval("20", 23103196, 23103238), 60, 484, 525, 2);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("523S970M"), true, new SimpleInterval("20", 23103633, 23104603), 60, 524, 1493, 3);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2, region3);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getChimericAlignmentsFromAlignmentRegions(contigSequence, alignmentRegionList, 50);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "contig-1");
        Assert.assertEquals(chimericAlignment.region1, region1);
        Assert.assertEquals(chimericAlignment.region2, region3);
        Assert.assertEquals(chimericAlignment.homology, "");
        Assert.assertEquals(chimericAlignment.insertedSequence, "TGAGAGTTGGCCCGAACACTGCTGGATTCCACTTCA");
        Assert.assertEquals(chimericAlignment.insertionMappings.size(), 1);
        Assert.assertEquals(chimericAlignment.insertionMappings.get(0), "1-contig-1:484-525:20,23103196,-,483S42M968S,60,2");
    }

    @Test
    public void testGetAssembledBreakpointFromAlignmentRegionsStrangeLeftBreakpoint() throws Exception {
        final byte[] contigSequence = LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, new SimpleInterval("chr19", 20138007, 20142231), 60, 1, contigSequence.length - 1986, 36);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, new SimpleInterval("chr19", 20152030, 20154634), 60, 3604, contigSequence.length, 36);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getChimericAlignmentsFromAlignmentRegions(contigSequence, alignmentRegionList, 50);
        Assert.assertEquals(assembledBreakpointsFromAlignmentRegions.size(), 1);
        final ChimericAlignment chimericAlignment = assembledBreakpointsFromAlignmentRegions.get(0);
        Assert.assertEquals(chimericAlignment.contigId, "702700");
        Assert.assertEquals(chimericAlignment.region1, region1);
        Assert.assertEquals(chimericAlignment.region2, region2);
        Assert.assertFalse(chimericAlignment.homology.isEmpty());

        final Tuple2<SimpleInterval, SimpleInterval> leftAndRightBreakpointsOnReferenceLeftAlignedForHomology = chimericAlignment.getLeftAndRightBreakpointsOnReferenceLeftAlignedForHomology();

        Assert.assertEquals(leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._1(), new SimpleInterval("chr19", 20138007, 20138007));
        Assert.assertEquals(leftAndRightBreakpointsOnReferenceLeftAlignedForHomology._2(), new SimpleInterval("chr19", 20152651, 20152651));
    }

    @Test
    public void testTreatAlignmentRegionAsInsertion() throws Exception {
        AlignmentRegion overlappingRegion1 = new AlignmentRegion("overlap", "22", TextCigarCodec.decode("47S154M"), false, new SimpleInterval("19", 48699881, 48700035), 60, 1, 154, 0);
        AlignmentRegion overlappingRegion2 = new AlignmentRegion("overlap", "22", TextCigarCodec.decode("116H85M"), true, new SimpleInterval("19", 48700584, 48700669), 60, 117, 201, 0);

        Assert.assertTrue(SVVariantCallerInternal.treatNextAlignmentRegionInPairAsInsertion(overlappingRegion1, overlappingRegion2, 50));
    }

    @Test
    public void testARtooSmall() {
        final byte[] contigSequence = LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, new SimpleInterval("chr19", 20138007, 20142231), 60, 1, contigSequence.length - 1986, 36);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, new SimpleInterval("chr19", 20152030, 20154634), 60, 3604, contigSequence.length, 36);
        Assert.assertFalse( SVVariantCallerInternal.alignmentRegionIsTooSmall(region1, region2, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );
        Assert.assertFalse( SVVariantCallerInternal.alignmentRegionIsTooSmall(region2, region1, SVConstants.DEFAULT_MIN_ALIGNMENT_LENGTH) );

        Assert.assertFalse( SVVariantCallerInternal.alignmentRegionIsTooSmall(region1, region2, 3000) );
        Assert.assertTrue( SVVariantCallerInternal.alignmentRegionIsTooSmall(region2, region1, 3000) );
    }

    @Test
    public void testGetHomology() {
        final byte[] contigSequence = "ATCGATCGAAAAGCTAGCTA".getBytes();

        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("12M8S"), true, new SimpleInterval("1", 1, 12), 60, 1, 12, 1);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("8H12M"), false, new SimpleInterval("1", 101, 112), 60, 9, 20, 1);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertEquals(SVVariantCallerInternal.getHomology(region2, region1, contigSequence), "AAAA");

        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("8M"), true, new SimpleInterval("1", 1, 12), 60, 1, 8, 1);            // dummy test data, almost guaranteed to be non-factual
        final AlignmentRegion region4 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("8M"), false, new SimpleInterval("1", 101, 112), 60, 13, 20, 1);    // dummy test data, almost guaranteed to be non-factual

        Assert.assertTrue(SVVariantCallerInternal.getHomology(region4, region3, contigSequence).isEmpty());
    }

    @Test
    public void testGetInsertedSequence() {
        final byte[] contigSequence = "GACGAACGATTTGACTTTAATATGAAATGTTTTATGTGGGCTATAAAATTATCCAAACTCGACACAGGACATTTTGAGCTTATTTCCAAATCATCTGGCCTTCATCTACCCACTGGAACTATTACTCTGCTGGGTCCTCATGGAAACATATCTTTCAGCCCTAACAATGAGACTACAGACATCTACGTCCCCAACACAACAGCTAAAAAGCAGTAGAATGTCAGAAAGGCTATCCACTTAGCCCTTGGCTGACAGGCCCCACTGAGCATCCTTTGCGAAGTCCATTTACTAGCTAATTCATAATTTACACAAGGCATTCAGACATAGCAGCTAAGATATAAAACATTTATCAACACAGGGACTAGTTTGTCATTTTAAAATAATTATGTTTAAGTAAGCCAATAAAGTCTATCTTCTCCAATTTACTTATTGAGCTTTATGAGGCAATTTAAGTCCCGATTTTGGGGGGTATGTATGAAAGGAGAGCATGGAAATGCCATTTGCTCCCTGAAGTTTTTATCTTTTTTTTTTTGAGATAGAGTCTTGTGTTTTCTGTGGAGTACATGAGTATGCATCAAAGCTAACAACGCCCACTGCCCTGTTAGTCAAATACCTTTGA".getBytes();
        final AlignmentRegion region1 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("532M87S"), true, new SimpleInterval("8", 118873207, 118873739), 60, 1, 532, 0);
        final AlignmentRegion region2 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("518S29M72S"), false, new SimpleInterval("1", 175705642, 175705671), 3, 519, 547, 0);
        final AlignmentRegion region3 = new AlignmentRegion("1", "contig-1", TextCigarCodec.decode("543S76M"), false, new SimpleInterval("1", 118875262, 118875338), 60, 544, 619, 0);

        Assert.assertTrue(SVVariantCallerInternal.getInsertedSequence(region1, region3, contigSequence).isEmpty());
        Assert.assertEquals(SVVariantCallerInternal.getInsertedSequence(region3, region1, contigSequence), "GAGATAGAGTC");

        Assert.assertTrue(SVVariantCallerInternal.getInsertedSequence(region1, region2, contigSequence).isEmpty() && SVVariantCallerInternal.getInsertedSequence(region2, region1, contigSequence).isEmpty());
    }

    // -----------------------------------------------------------------------------------------------
    // Inversion specific
    // -----------------------------------------------------------------------------------------------

    @Test
    public void testProduceAlleles_Inversion() throws IOException {
        final PipelineOptions options = null;
        final ReferenceMultiSource referenceMultiSource = new ReferenceMultiSource(options, twoBitRefURL, ReferenceWindowFunctions.IDENTITY_FUNCTION);
        final List<Allele> alleles = SVVariantCallerInternal.produceAlleles(referenceMultiSource, "20", 1000000, 2000000);
        Assert.assertEquals(alleles.size(), 2);
        Assert.assertTrue(alleles.get(0).isReference() && alleles.get(1).isNonReference());
        Assert.assertTrue(alleles.get(1).isSymbolic());
    }

    @Test
    public void testProduceVariantId_Inversion() {
        final byte[] contigSequence = LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, new SimpleInterval("chr19", 20138007, 20142231), 60, 1, contigSequence.length - 1986, 36);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, new SimpleInterval("chr19", 20152030, 20154634), 60, 3604, contigSequence.length, 36);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getChimericAlignmentsFromAlignmentRegions(contigSequence, alignmentRegionList, 50);

        final BreakpointAlleleInversion inversionAllele = new BreakpointAlleleInversion(assembledBreakpointsFromAlignmentRegions.get(0));
        final String id = SVVariantCallerInternal.produceVariantId(inversionAllele);
        Assert.assertFalse(id.isEmpty());
        final String[] fields = id.split(SVConstants.VARIANT_ID_FIELD_SEPARATER);
        Assert.assertEquals(fields.length, 4);
        Assert.assertEquals(BreakpointAlleleInversion.InversionType.valueOf(fields[0]), BreakpointAlleleInversion.InversionType.INV_3_TO_5);
        Assert.assertTrue(fields[1].endsWith("19"));
        Assert.assertEquals(fields[2], "20138007");
        Assert.assertEquals(fields[3], "20152651");
    }

    @Test
    public void testUpdateAttributes_Inversion() {
        // not an exhaustive test on all attributes
        // MAPPING_QUALITIES, ALIGNMENT_LENGTH, INSERTED_SEQUENCE, INSERTED_SEQUENCE_MAPPINGS, HOMOLOGY, TYPE

        final byte[] contigSequence = LONG_CONTIG1.getBytes();
        AlignmentRegion region1 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("1986S236M2D1572M1I798M5D730M1I347M4I535M"), false, new SimpleInterval("chr19", 20138007, 20142231), 60, 1, contigSequence.length - 1986, 36);
        AlignmentRegion region2 = new AlignmentRegion("702700", "702700", TextCigarCodec.decode("3603H24M1I611M1I1970M"), true, new SimpleInterval("chr19", 20152030, 20154634), 60, 3604, contigSequence.length, 36);
        final List<AlignmentRegion> alignmentRegionList = Arrays.asList(region1, region2);
        final List<ChimericAlignment> assembledBreakpointsFromAlignmentRegions = SVVariantCallerInternal.getChimericAlignmentsFromAlignmentRegions(contigSequence, alignmentRegionList, 50);

        final Iterable<Tuple3<String, String, ChimericAlignment>> alignments = Arrays.asList(new Tuple3<>("702700", "702700", assembledBreakpointsFromAlignmentRegions.get(0)));

        final BreakpointAllele inversionAllele = new BreakpointAlleleInversion(assembledBreakpointsFromAlignmentRegions.get(0));

        VariantContextBuilder vcBuilder = new VariantContextBuilder().chr("chr19").start(20138007).stop(20152651).alleles(Arrays.asList(Allele.create("A", true), Allele.create("<INV>"))); // dummy test data, almost guaranteed to be non-factual
        final VariantContext vc = SVVariantCallerInternal.updateAttributes(vcBuilder, inversionAllele, alignments, 20138007, 20152651).make();

        String testString = vc.getAttributeAsString(GATKSVVCFHeaderLines.MAPPING_QUALITIES, "NONSENSE");
        Assert.assertFalse(testString.isEmpty());
        Assert.assertEquals(testString.split(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR), new String[]{"60"});

        String[] ss = vc.getAttributeAsString(GATKSVVCFHeaderLines.ALIGN_LENGTHS, "NONSENSE").split(GATKSVVCFHeaderLines.FORMAT_FIELD_SEPARATOR);
        Assert.assertEquals(ss.length, 1);
        Assert.assertEquals((int)Integer.valueOf(ss[0]), Math.min(region1.referenceInterval.size(), region2.referenceInterval.size()) - SVVariantCallerUtils.overlapOnContig(region1, region2));

        Assert.assertTrue(vc.getAttributeAsString(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, "").isEmpty());
        Assert.assertTrue(vc.getAttributeAsString(GATKSVVCFHeaderLines.INSERTED_SEQUENCE_MAPPINGS, "").isEmpty());

        Assert.assertFalse(vc.getAttributeAsString(GATKSVVCFHeaderLines.HOMOLOGY, "").isEmpty()); // simple test, as more complete test is covered in {@link testGetHomology()}

        Assert.assertEquals(vc.getAttributeAsString(GATKSVVCFHeaderLines.SVTYPE, "NONSENSE"), GATKSVVCFHeaderLines.SVTYPES.INV.name());
    }

    // -----------------------------------------------------------------------------------------------
    // Data providers
    // -----------------------------------------------------------------------------------------------

    private static String twoBitRefURL = publicTestDir + "large/human_g1k_v37.20.21.2bit";

    static final String LONG_CONTIG1 = "TTTTTTTTTTTTTTTCTGAGACCGAATCTCGCTCTGTCACCCAGGCTGGAGTGCAGTGGCACGATCTTGGCTTACTGCAAGCTCTGCCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCCCACCCCCCCACCCCCCCAGGTAGCTG" +
            "GGACTACAGGTGTCTGCCACCACACCTGGCTAAATTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAAGATGGTTTCGCTCTCCTGACCTCGCGATCCGCCCACCTCGGCCTCTCAAAGTGCTGGGATTACAGGCCTGAGCCACTGCGCCC" +
            "TGCCTGACCCCTTCTTTTAAAACAAATCTTTTGTCTTTGTCTTCATTTCTGCATTCGTCCCCTTCGTTCAATCCTGTAGGAATTGACAGTGATATTGGGGAGTTCCCATCCCTGGATTTGGGATTTCCTCGAGTTTCCAGCCCTGTCCTTGTGGCCAAAAAGT" +
            "GGCACAGTTTGGTTTTTGAATCCTGTTCGTGACGCCAAAAATATTCTGTGTGGGAAACATTCAAGGGGAGAAGAAAAGACACACACACAATACCTTTAAGGGTAAATAAGCTTTGTGCCACGTAAATGGAAATGCAGATATAATTAGCAAATTATATAATAAG" +
            "CAAATCAATGTAATAACCAGATTGATATAATAAGCATGTTGATATAATAAGCAAATTGCAATGGGAAGAGGAGAAAGGAAAAGAGATATATATATTTACACTCACCAGACTATGGAGGATTCACCACCAGACTGGGAAGCAACAGCTTGGGCTCCAGAGTCAG" +
            "CTACTCATTCCTGCACAGATGAGGAGGGTCTAATGAAGCTTCAGCACAATCTGATACCCTAGCTCTTTTGTAACGAGTTGTTTGGCATAAGGCCCAGTCATGAGGGCCATTTGCAACTGGGGTCAAGGAACACAAAACTGTCAACTTGTTTTTGCGATTGTCT" +
            "ATTGTTTTTCAACAACTAATATATAGAAATAGATTGAAATAGAGATTTCTCTGAAACAGCGCTGAATGAATGCCTCAAGGGCCTCACACAACCTGTTTGAGAACTTGGTGACTACTGAGTTTGTCCACGTTCAATTAAGTTCAAATTTAGTATTTAACTTTTC" +
            "CTCCACAAATTACCCAGTCTTGATGAATCAGCTCTGTCTAGACATAGGGCAAGGTGAACCCCTTGGGCAGTTACACAACCTCCGCCTTCTGGGTTTAAGCAATTCTCCTGCCTCAGCCTCCGGACTAGCTGGGTCTACAGGTGTGCAGCACCACACCCAGCTA" +
            "GTTATTTGTACTTTTAGTAGAAATGGGGTTTCACCATGTTGGCCAGGCTGGTCTTGAACTCCTGACCTCAAGTGATCCACCCACCTTGGCCTCCCAAAGTGTTGCGATTACAGGCACTTGCCAGTGAACCTGGCCCTAAATGACTTCTTTCTATCTCCTATAA" +
            "CAATTTGAAATTACTTAAAGGTGGTTTCAAATTGAAAAAATAAAAAGAATTTGGATAAAAATAAAATGTAAACAGTTTTTAAAAATTACAAGAGATTACAAAATATACATGTAAAACCGGAGTGGTCAAAAATGACAAATTTGATTTATTTATAAGGTTTATT" +
            "AAAATTAGCTTTAGTATTGATAATACACTATTACCAAAGTAAAAGCTGATTTTCTCTTGAAAAAAATTTTATGTATTATTAATATGACAGCAAAATACTTCTGTTCACCTTTTGAATATATTCAAAAAGAGAGAGAGTAAAAAAAAAGATAGAATTTTCCCAT" +
            "GATCTGGGTTGGGCCTGGCTTAGCTCAGGGAGGAAGCCCTGCCTGAAAAACGCTGCAGCTTAGGCTGTGACTTTTTCTTCACTCAGTCCTGATCACATTTTCTGTTACTCAGGGCCTGAGGGGGCGGGGGCCTTAAGCATTATCCAATCAGAAACGCTGGGCT" +
            "GACGCCCCGTCCGGGAGGGAGGTGGGGGGGGTCAGCCCCCCGCCCAGCCAGCCGCCCCGTCCGGGAGGTGGGGGGTGCCTCTGCCCGGCCGCCCCTACTGGGAAGTGAGGAGCCCCTCTGCCCGGCCACCACCCCGTCTGGGAGGTGTACCCAACAGCTCATT" +
            "GAGAACGGGCCATGATGACAATGGCAGTTTTGTGGAATAGAAACGGGGGAAAGGTGGGGAAAAGATTGAGAAATCGGATGGTTGCTGTGTCTGTGTAGAAAGAGGTAGACATGGGAGACTTCATTTTGTTCTGTACTAAGAAAAATTCTTCTGCCTTGGGATC" +
            "CTGTTGACCTATGACCTTACCCCCAACCCTGTGCTCTCTGAAACATGTGCTGTGTCCACTCAGGGTTAAATGGATTAAGGGCGGTGCAAGATGTGCTTTGTTAAACAGATGCTTGAAGGCAGCATGCTCGTTAAGAGTCATCACCACTCCCTAATCTCAAGTA" +
            "CCCGGGGACACAAACACTGCGGAAGGCCGCAGGGTCCTCTGCCTAGGAAAACCAGAGACCTTTGTTCACTTGTTTATCTGCTGACCTTCCCTCCACTGTTGTCCTATGACCCTGCCAAATCCCCCTCTGCGAGAAACACCCAAGAATGATCAACTAAAAAAAA" +
            "AAAAGAAAAAAGAAAAAAGAAAAAATACACCCTGGCAGAAAACTTCTCGAATGTAAGTAAAGTGGGAGTTTTCAACTTGCCCACCTCCTTACAAAAACATGCAATAAGAATGTGCTAGATATAAAGCTTACAAATGTGAGAAATACGGGAAGTTTTTCAATTT" +
            "TAACAAATACTTTCAAAGTCATGTGAAAATTCTCACTAGAAAGAAATCCAATAAAGGCAAGTTGAGGCATAAAAAAAAAAAAAAGAAAAAGCTCAGAATACATATAAGATTGACTAAGTCAGCCAGAAAATAATCCCCTAAAAGAAATTTCTCTCTAAACACC" +
            "CAATGTGCACAGCTACTCTCAGCATGAGAAACATGAGCTTTATGAAGAAAGGTGGCAGATTTTCAGAAGAATTTTATAAAAGTTTCTTTTCCATCTCTGCTGCTCTCTCATCTCCTAGCCATTGAATGGGGGTTCTATATTGAAATACATCTGACAACTTCCA" +
            "ACAACACTTTTTGATCAAGAAATAGAATTTGACTATGTTCGTATAGTGGAATATATTAGAACTTGTAACACAGCTAACTGAATAGCTATTATGGTGTTTGGGTGGCCACATCACCTGTCTTTATTTGTCCGGTAATAGCAGCATTCCAATTTAAAGAAATAAA" +
            "AGATACCAAAATTGTGTTTACTTTTAATTATTCCTATTGAATAAAGTAATAAGCATGTCAGACTGATATCTATCATAACAATAAATTTTGTTTGGATATTATATTAGATATAAATATTTAAGTATGAATAATTTTAATGAACTAGTCATAATGTATGTAGCAT" +
            "TTTTAAAAATTGTAACTATACTTCAGTTAAAACACTTTATATTTCAAAAGCATAAATAACAATATTAAAATAACAATTTAGGTGATTCATTCAAAGTAAGTATTGGGGCTTTATGTTCATACTATTGTAGAAAATACTGCTTATGGCTCATGCCTGTAATCCC" +
            "AGCACATTGGGAGGCTGAGGTGGGTAGATCACCTGAGGTCAGGAGTTCCTGATCCCATCTTTACTAAAAATACAAAACTTACCCAGGGTGGTTGTGCACACTTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATTGCTTGAACAAGGGAGGAAATGG" +
            "TTGCAGTGAGCCATGATCATGCCACTGAACCCCAGCCTGGGCAAGAGAGTGAGACTGTCTCAAAAAAAAAAAAAACTGTTTAATTTTTATGAATGCAGGTTTTCTGCAAACACTACACATAACTATGCTAATTGTTCTGAAGTAATAAATAGAAAGCAAGGCA" +
            "CAACTACAGACTCCACTGTTCAGTTTATGCACTGAACTGTTCTTGCTTTTGCAGTGTAAGTATTTCTGCCTGCAAATACTGGATAATTACCTTGGATCATCAGATTTCTATCAAAGGAATTTAGTATCTTTTAGTCTTTATCATTTTGTATTGCTAAATTTAT" +
            "CTGTGTGTTAAGCTTCTGTGTGCTCTTAAAATGAGGTTTTATCTAAACAAACCTGTGTCTACTTTAAAAGACTAAACATGAAAAAACTAAACTTTTCAGAACCAAAAACAAAGCAATAAATCTGAAGTACTAGATAGTCTGGAGTGAGATTTATTTAGCTTTT" +
            "TTTTTTTTTTGAGATGGAGTCTCGCTCTGTCACCGAGGCTGGAGTGCAGTGGCACGAACTCGGCTCACTGCAAAAGCTCTGCCTCCCAGCTTCATGCCATTCTCCTACCTCAGCCTCCCAAGTAGCTGGGATTACAGGCAACTGCCACCACGCCCAGCTAATT" +
            "TTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTTGATCTCCTGACCTCATGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCATGCAAAACCGCGCCCTGCCCTTATTTAGCTCTTAATATGATTTACATATAT" +
            "TTCAAAAAAGCAGAGAAAAATATCTACATATAATCTAAATCCCTTAAGAGAGAATAGCAAAAAAATTTTTGGAACTTTTTAAAGAGTTTTGAACTCTTGGACATCTGAATTTTGCACACTATATGCACTTGAAAGAATGTTAATGGGGAAAAAGCAGAAGAGA" +
            "AAAAGACGTTATAAAAAAATCCATGAGTGCACAACACCTATGAAACAGAATAGAGAGCCCAGAAATAATGCCTTTCACCTGCAACCATCAGATTTCTGACAAAGCTGACAAGAGGAATGTGGGAAGAATTCTCTCTTTCATAAATGGTGCTGGAATAACTATC" +
            "TACCACTATGTAGAAGACTGAAGTGGACCCCTTCATTACACCATATAAAAAAATCAACTGAAGATAAATTAAGGACTTAAATGTAAAACTTAAAATTATAAGAAACCCTGCAAGATAACCTAGGAAATAGCATTCTAGACACAGAAACAGGTAAAGACTTCAT" +
            "GATGAAGCTACCAAAAGCAACTGCAACAGAAGTAAATTGACAAATGGGATGTATTTAAACTTAAGAGCTTCTTCACAGCAAAGGAAACTATCAACAGAGTAAACAGACAAACTAGAGAATAAAAGAATATATTTGTAAATTTTGCCTCTGAAAAAGGTCTAAT" +
            "ATACAGAATTTATTAGGAACTTAAACAAGTTTACAAGAAAAAAACACACTCATTAAAAAGTATGCAAAAAACATGAACAGATGCCTTTCAATAGAAGATACACATAGGGCTAACAAGCATATGAAAAAAAAAATGCTCATTGCTAATCATTAGAGAAATTAGA" +
            "AGAGAAATCACAGGCTGGGTGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCAAGGCAGGCAGATTACAAGGTCAGGAGATCAAGATCATCCTGGCTAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAAATCAGCCAGGTGTGGCAGTGT" +
            "GCACCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATTGCTTGAATCTGGTAGGCAGAGGTTGCAGTGAGCTGAGATCACACCACTGCACTCCTGCCTGGGCAACAGAGCAAGACTCCGTCTCAAACACACACACACAGACACACACACACACACAC" +
            "ACACACACACACACACACACGCAGAGAAACCACAATGAGATACCACCTCATACCAGTCAGAATGGCTATTTTTAAAAAGTCAAAAGATAACAGATGCTGACAAAGTTGCAGAGAAAAGGGAATGCTTATACTCCTCTGGTGGGAGTGTAAATTATTTCAACAA" +
            "CTGTAAAAATCAGTGTGACAATTCCTCACAGAATGAAAAACAGAATTATCATTCGACTGGGAAACCTCATAATTGAGTATATACCCAAAGAAATATGAAATATTATAAAGACACATCCACATGCATGTTCACTGCAGCACTATTCACAATAGCAAAGACACGG" +
            "ACAGACTAAATGCCTATCAATGGCAGACTGGATCAAGAAAATATGGTATGGTCAGATGCGGTGGCTCATGCCTGTAATTCCAGCCCTTTGGGAGGCTGAGGCAGGTGGATTGCCTGAGCTTAGAAGTTTGAGACCACTCTGGGCAACATGGCAAAATTTTGTC" +
            "TCCACAGAAGATACAAAAAAAAAAAAAAAAAA";
}