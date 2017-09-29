package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.base.Strings;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SWPairwiseAlignment;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link HaplotypeCallerGenotypingEngine}.
 */
public final class HaplotypeCallerGenotypingEngineUnitTest extends BaseTest {

    private class BasicGenotypingTestProvider extends TestDataProvider {
        byte[] ref;
        byte[] hap;
        Map<Integer,Byte> expected;

        public BasicGenotypingTestProvider(String refString, String hapString, Map<Integer, Byte> expected) {
            super(BasicGenotypingTestProvider.class, String.format("Haplotype to VCF test: ref = %s, alignment = %s", refString,hapString));
            ref = refString.getBytes();
            hap = hapString.getBytes();
            this.expected = expected;
        }
        
        public Map<Integer,VariantContext> calcAlignment() {
            final SWPairwiseAlignment alignment = new SWPairwiseAlignment(ref, hap, new SWPairwiseAlignment.Parameters(3,-1,-4, -1));
            final Haplotype h = new Haplotype(hap, false, alignment.getAlignmentStart2wrt1(), alignment.getCigar());
            return new EventMap(h, ref, new SimpleInterval("4", 1, 1 + ref.length), "name");
        }

        public String toString() {
            return "REF:" + new String(ref) + ",ALT:" + new String(hap);
        }
    }

    @DataProvider(name = "BasicGenotypingTestProvider")
    public Object[][] makeBasicGenotypingTests() {

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(2 + contextSize, (byte)'M');
            map.put(21 + contextSize, (byte)'M');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG", "ATCTCGCATCGCGAGCATCGCCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1, (byte)'M');
            map.put(20, (byte)'D');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider("AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(2 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            map.put(30 + contextSize, (byte)'D');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "ACCTCGCATCGCGAGCATCGTTACTAGCCGATG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            map.put(28 + contextSize, (byte)'M');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCCATAG", map);
        }

        return BasicGenotypingTestProvider.getTests(BasicGenotypingTestProvider.class);
    }
    
    @Test(dataProvider = "BasicGenotypingTestProvider", enabled = true)
    public void testHaplotypeToVCF(BasicGenotypingTestProvider cfg) {
        Map<Integer,VariantContext> calculatedMap = cfg.calcAlignment();
        Map<Integer,Byte> expectedMap = cfg.expected;
        logger.warn(String.format("Test: %s", cfg.toString()));
        if(!compareVCMaps(calculatedMap, expectedMap)) {
            logger.warn("calc map = " + calculatedMap);
            logger.warn("expected map = " + expectedMap);
        }
        Assert.assertTrue(compareVCMaps(calculatedMap, expectedMap),"" + cfg);
    }

    @Test(dataProvider="AddMiscellaneousDataProvider", enabled=false)
    public void testAddMiscellaneousAllele(final String readBases, final int readOffset,
                                           final String ref, final int refOffset,
                                           final String referenceAllele, final String[] alternatives, final double[] likelihoods, final double[] expected) {
        final byte baseQual = (byte)30;

        final byte[] baseQuals = Utils.dupBytes(baseQual, readBases.length());
        final GATKRead read = ArtificialReadUtils.createArtificialRead(readBases.getBytes(), baseQuals, readBases.length() + "M");
        final Locatable loc = new SimpleInterval("20",refOffset,refOffset);
        final ReadPileup pileup = new ReadPileup(loc, Collections.singletonList(read),readOffset);
        final VariantContextBuilder vcb = new VariantContextBuilder();
        final GenotypeBuilder gb = new GenotypeBuilder();
        final List<String> alleleStrings = new ArrayList<>( 1 + alternatives.length);
        alleleStrings.add(referenceAllele);
        alleleStrings.addAll(Arrays.asList(alternatives));

        gb.AD(new int[] { 1 });
        gb.DP(1);
        gb.PL(likelihoods);

        vcb.alleles(alleleStrings);
        vcb.loc("20",refOffset,refOffset + referenceAllele.length() -1);

        vcb.genotypes(gb.make());

        final VariantContext vc = vcb.make();

        final VariantContext updatedVc = null; // GenotypingEngine.addMiscellaneousAllele(vc,pileup,ref.getBytes(),0);
        final GenotypeLikelihoods updatedLikelihoods = updatedVc.getGenotype(0).getLikelihoods();
        Assert.assertEquals(updatedLikelihoods.getAsVector().length, expected.length);
        final double[] updatedLikelihoodsArray = updatedVc.getGenotype(0).getLikelihoods().getAsVector();
        for (int i = 0; i < updatedLikelihoodsArray.length; i++) {
            Assert.assertEquals(updatedLikelihoodsArray[i],expected[i],0.0001);
        }
        Allele altAllele = null;
        for (final Allele allele : updatedVc.getAlleles())
            if (allele.isSymbolic() && allele.getBaseString().equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME))
                altAllele = allele;
        Assert.assertNotNull(altAllele);
    }

    @DataProvider(name="AddMiscellaneousDataProvider")
    public Iterator<Object[]> addMiscellaneousAlleleDataProvider() {
        return Arrays.asList(ADD_MISCELLANEOUS_ALLELE_DATA).iterator();
    }

    private static final double MATCH_LnLK = QualityUtils.qualToProbLog10((byte)30);
    private static final double MISS_LnLK = QualityUtils.qualToErrorProbLog10((byte)30);

    private static final Object[][] ADD_MISCELLANEOUS_ALLELE_DATA = new Object[][] {
            new Object[] {"ACTG", 0,"ACTGTGAGTATTCC",0,"A",new String[]{}, new double[] {MATCH_LnLK * MATCH_LnLK}, 6,
                    new double[] {MATCH_LnLK * MATCH_LnLK,MATCH_LnLK * MISS_LnLK, MISS_LnLK * MISS_LnLK}}
    };

    /**
     * Private function to compare Map of VCs, it only checks the types and start locations of the VariantContext
     */
    private boolean compareVCMaps(Map<Integer, VariantContext> calc, Map<Integer, Byte> expected) {
        if( !calc.keySet().equals(expected.keySet()) ) { return false; } // sanity check
        for( Integer loc : expected.keySet() ) {
            Byte type = expected.get(loc);
            switch( type ) {
                case 'I':
                    if( !calc.get(loc).isSimpleInsertion() ) { return false; }
                    break;
                case 'D':
                    if( !calc.get(loc).isSimpleDeletion() ) { return false; }
                    break;
                case 'M':
                    if( !(calc.get(loc).isMNP() || calc.get(loc).isSNP()) ) { return false; }
                    break;
                default:
                    return false;
            }
        }
        return true;
    }


    @DataProvider(name = "CreateHaplotypeMappingProvider")
    public Object[][] makeCreateHaplotypeMappingData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Set<Haplotype> haplotypes = new HashSet<>();
        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);

        final Haplotype AtoC1 = new Haplotype("AACAA".getBytes());
        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altC)).make();
        AtoC1.setEventMap(new EventMap(Arrays.asList(vc1)));
        AtoC1.getEventMap().put(3, vc1);
        haplotypes.add(AtoC1);

        final Haplotype AtoC2 = new Haplotype("AAACA".getBytes());
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(4).stop(4).alleles(Arrays.asList(ref, altT)).make();
        AtoC2.setEventMap(new EventMap(Arrays.asList(vc2)));
        AtoC2.getEventMap().put(4, vc2);
        haplotypes.add(AtoC2);

        tests.add(new Object[]{vc1, haplotypes, AtoC1});
        tests.add(new Object[]{vc2, haplotypes, AtoC2});
        tests.add(new Object[]{new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altT)).make(), haplotypes, null});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="CreateHaplotypeMappingProvider")
    public void testCreateHaplotypeMapping(final VariantContext vc, final Set<Haplotype> haplotypes, final Haplotype expected) {
        final Map<VariantContext, Set<Haplotype>> mapping = HaplotypeCallerGenotypingEngine.constructHaplotypeMapping(Arrays.asList(vc), haplotypes);
        final Set<Haplotype> actual = mapping.get(vc);
        if ( expected == null )
            Assert.assertTrue(actual.isEmpty(), actual.toString());
        else {
            Assert.assertEquals(actual.size(), 1);
            Assert.assertEquals(actual.iterator().next(), expected);
        }
    }

    @DataProvider(name = "ConstructPhaseSetMappingProvider")
    public Object[][] makeConstructPhaseSetMappingData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);

        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(2).stop(2).alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc3 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altT)).make();
        final VariantContext vc4 = new VariantContextBuilder().chr("20").start(4).stop(4).alleles(Arrays.asList(ref, altC)).make();
        final List<VariantContext> calls = Arrays.asList(vc2, vc3, vc4);

        final Haplotype pos1 = new Haplotype("CAAAA".getBytes());
        pos1.setEventMap(new EventMap(Arrays.asList(vc1)));
        pos1.getEventMap().put(1, vc1);
        final Haplotype pos2 = new Haplotype("ACAAA".getBytes());
        pos2.setEventMap(new EventMap(Arrays.asList(vc2)));
        pos2.getEventMap().put(2, vc2);
        final Haplotype pos3 = new Haplotype("AACAA".getBytes());
        pos3.setEventMap(new EventMap(Arrays.asList(vc3)));
        pos3.getEventMap().put(3, vc3);
        final Haplotype pos4 = new Haplotype("AAACA".getBytes());
        pos4.setEventMap(new EventMap(Arrays.asList(vc4)));
        pos4.getEventMap().put(4, vc4);
        final Haplotype pos24 = new Haplotype("ACACA".getBytes());
        pos24.setEventMap(new EventMap(Arrays.asList(vc2, vc4)));
        pos24.getEventMap().put(2, vc2);
        pos24.getEventMap().put(4, vc4);
        final Haplotype pos34 = new Haplotype("AACCA".getBytes());
        pos34.setEventMap(new EventMap(Arrays.asList(vc3, vc4)));
        pos34.getEventMap().put(3, vc3);
        pos34.getEventMap().put(4, vc4);
        final Haplotype pos234 = new Haplotype("ACCCA".getBytes());
        pos234.setEventMap(new EventMap(Arrays.asList(vc2, vc3, vc4)));
        pos234.getEventMap().put(2, vc2);
        pos234.getEventMap().put(3, vc3);
        pos234.getEventMap().put(4, vc4);

        final Map<VariantContext, Set<Haplotype>> haplotypeMap = new HashMap<>();

        // test no phased variants #1
        final Set<Haplotype> haplotypes2 = new HashSet<>();
        haplotypes2.add(pos2);
        haplotypeMap.put(vc2, haplotypes2);
        tests.add(new Object[]{Arrays.asList(vc2), new HashMap<>(haplotypeMap), 2, 0, 0, 0, 0});

        // test no phased variants #2
        final Set<Haplotype> haplotypes3 = new HashSet<>();
        haplotypes3.add(pos3);
        haplotypeMap.put(vc3, haplotypes3);
        tests.add(new Object[]{Arrays.asList(vc2, vc3), new HashMap<>(haplotypeMap), 3, 0, 0, 0, 0});

        // test opposite phase
        tests.add(new Object[]{Arrays.asList(vc2, vc3), new HashMap<>(haplotypeMap), 2, 2, 1, 1, 1});

        // test no phased variants #3
        final Set<Haplotype> haplotypes4 = new HashSet<>();
        haplotypes4.add(pos4);
        haplotypeMap.put(vc4, haplotypes4);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 3, 0, 0, 0, 0});

        // test mixture
        final Set<Haplotype> haplotypes24 = new HashSet<>();
        haplotypes24.add(pos24);
        haplotypeMap.put(vc2, haplotypes24);
        haplotypeMap.put(vc4, haplotypes24);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 2, 1});

        // test 2 hets
        haplotypeMap.remove(vc3);
        tests.add(new Object[]{Arrays.asList(vc2, vc4), new HashMap<>(haplotypeMap), 1, 2, 1, 2, 0});

        // test 2 with opposite phase
        final Set<Haplotype> haplotypes1 = new HashSet<>();
        haplotypes1.add(pos1);
        haplotypeMap.put(vc1, haplotypes1);
        tests.add(new Object[]{Arrays.asList(vc1, vc2, vc4), new HashMap<>(haplotypeMap), 2, 3, 1, 1, 2});

        // test homs around a het
        final Set<Haplotype> haplotypes2hom = new HashSet<>();
        haplotypes2hom.add(pos24);
        haplotypes2hom.add(pos234);
        final Set<Haplotype> haplotypes4hom = new HashSet<>();
        haplotypes4hom.add(pos24);
        haplotypes4hom.add(pos234);
        final Set<Haplotype> haplotypes3het = new HashSet<>();
        haplotypes3het.add(pos234);
        haplotypeMap.put(vc2, haplotypes2hom);
        haplotypeMap.put(vc3, haplotypes3het);
        haplotypeMap.put(vc4, haplotypes4hom);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 3, 0});

        // test hets around a hom
        final Set<Haplotype> haplotypes2het = new HashSet<>();
        haplotypes2het.add(pos234);
        final Set<Haplotype> haplotypes4het = new HashSet<>();
        haplotypes4het.add(pos234);
        final Set<Haplotype> haplotypes3hom = new HashSet<>();
        haplotypes3hom.add(pos3);
        haplotypes3hom.add(pos234);
        haplotypeMap.put(vc2, haplotypes2het);
        haplotypeMap.put(vc3, haplotypes3hom);
        haplotypeMap.put(vc4, haplotypes4het);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 2, 3, 1, 3, 0});

        // test no phased variants around a hom
        final Set<Haplotype> haplotypes2incomplete = new HashSet<>();
        haplotypes2incomplete.add(pos24);
        final Set<Haplotype> haplotypes3incomplete = new HashSet<>();
        haplotypes3incomplete.add(pos34);
        final Set<Haplotype> haplotypes4complete = new HashSet<>();
        haplotypes4complete.add(pos24);
        haplotypes4complete.add(pos34);
        haplotypes4complete.add(pos234);
        haplotypeMap.put(vc2, haplotypes2incomplete);
        haplotypeMap.put(vc3, haplotypes3incomplete);
        haplotypeMap.put(vc4, haplotypes4complete);
        tests.add(new Object[]{calls, new HashMap<>(haplotypeMap), 0, 0, 0, 0, 0});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="ConstructPhaseSetMappingProvider")
    public void testConstructPhaseSetMapping(final List<VariantContext> calls,
                                             final Map<VariantContext, Set<Haplotype>> haplotypeMap,
                                             final int totalHaplotypes,
                                             final int expectedMapSize,
                                             final int expectedNumGroups,
                                             final int expectedNum01,
                                             final int expectedNum10) {
        final Map<VariantContext, Pair<Integer, String>> actualPhaseSetMapping = new HashMap<>();
        final int actualNumGroups = HaplotypeCallerGenotypingEngine.constructPhaseSetMapping(calls, haplotypeMap, totalHaplotypes, actualPhaseSetMapping);
        Assert.assertEquals(actualNumGroups, expectedNumGroups);
        Assert.assertEquals(actualPhaseSetMapping.size(), expectedMapSize);

        int num01 = 0, num10 = 0;
        for ( final Pair<Integer, String> phase : actualPhaseSetMapping.values() ) {
            if ( phase.getRight().equals("0|1") )
                num01++;
            else if ( phase.getRight().equals("1|0") )
                num10++;
        }
        Assert.assertEquals(num01, expectedNum01);
        Assert.assertEquals(num10, expectedNum10);
    }

    @DataProvider(name = "ConstructPhaseGroupsProvider")
    public Object[][] makeConstructPhaseGroupsData() {
        List<Object[]> tests = new ArrayList<>();

        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);

        final Genotype g1 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc1 = new VariantContextBuilder().chr("20").start(1).stop(1).alleles(Arrays.asList(ref, altC)).genotypes(g1).make();
        final Genotype g2 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc2 = new VariantContextBuilder().chr("20").start(2).stop(2).alleles(Arrays.asList(ref, altC)).genotypes(g2).make();
        final Genotype g3 = new GenotypeBuilder().alleles(Arrays.asList(ref, altC)).make();
        final VariantContext vc3 = new VariantContextBuilder().chr("20").start(3).stop(3).alleles(Arrays.asList(ref, altC)).genotypes(g3).make();
        final List<VariantContext> calls = Arrays.asList(vc1, vc2, vc3);

        // test no phased variants, empty map
        final Map<VariantContext, Pair<Integer, String>> nonePhased1 = new HashMap<>();
        tests.add(new Object[]{calls, nonePhased1, 0, 0, 0});

        // test no phased variants, full map, exception expected
        final Map<VariantContext, Pair<Integer, String>> nonePhased2 = new HashMap<>();
        nonePhased2.put(vc1, Pair.of(0, "0/1"));
        nonePhased2.put(vc2, Pair.of(1, "0/1"));
        nonePhased2.put(vc3, Pair.of(2, "0/1"));
        tests.add(new Object[]{calls, nonePhased2, 3, -1, -1});

        // test 2 phased variants
        final Map<VariantContext, Pair<Integer, String>> twoPhased = new HashMap<>();
        twoPhased.put(vc1, Pair.of(0, "0/1"));
        twoPhased.put(vc2, Pair.of(0, "0/1"));
        tests.add(new Object[]{calls, twoPhased, 1, 1, 2});

        // test all phased variants
        final Map<VariantContext, Pair<Integer, String>> allPhased = new HashMap<>();
        allPhased.put(vc1, Pair.of(0, "0/1"));
        allPhased.put(vc2, Pair.of(0, "0/1"));
        allPhased.put(vc3, Pair.of(0, "0/1"));
        tests.add(new Object[]{calls, allPhased, 1, 1, 3});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider="ConstructPhaseGroupsProvider")
    public void testConstructPhaseGroups(final List<VariantContext> calls,
                                         final Map<VariantContext, Pair<Integer, String>> phaseMap,
                                         final int endIndex,
                                         final int expectedNumGroups,
                                         final int expectedGroupSize) {
        final List<VariantContext> actualPhasedCalls;
        try {
            actualPhasedCalls = HaplotypeCallerGenotypingEngine.constructPhaseGroups(calls, phaseMap, endIndex);
        } catch (IllegalStateException e) {
            Assert.assertEquals(-1, expectedNumGroups);
            return;
        }

        final Set<String> uniqueGroups = new HashSet<>();
        int counter = 0;
        for ( final VariantContext call : actualPhasedCalls ) {
            for ( final Genotype g : call.getGenotypes() ) {
                if ( g.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY) ) {
                    uniqueGroups.add(g.getExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY).toString());
                    counter++;
                }
            }
        }

        Assert.assertEquals(uniqueGroups.size(), expectedNumGroups);
        Assert.assertEquals(counter, expectedGroupSize);
    }

    @Test
    public void testReduceNumberOfAlternativeAllelesBasedOnHaplotypesScores(){

        // first have a list of alleles, one ref, several alt
        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);
        final Allele altT2 = Allele.create("TT", false);
        final Allele altG = Allele.create("G", false);

        // then create several haplotypes, assign ad-hoc scores
        final Haplotype hapRef = new Haplotype("AAAAA".getBytes());
        hapRef.setScore(Double.MAX_VALUE);

        // test case when both same best score and second best score are the same
        final Haplotype hapT = new Haplotype("TAAAA".getBytes());
        hapT.setScore(-2.0);
        final Haplotype hapTAnother = new Haplotype("TAAAT".getBytes());
        hapTAnother.setScore(-3.0);
        final Haplotype hapT2 = new Haplotype("TTAAA".getBytes());
        hapT2.setScore(-2.0);
        final Haplotype hapT2Another = new Haplotype("TTAAT".getBytes());
        hapT2Another.setScore(-3.0);

        final Haplotype hapC = new Haplotype("CAAAA".getBytes());
        hapC.setScore(-3.0);

        // for case when there's tie in highest haplotype score
        final Haplotype hapG = new Haplotype("GAAAA".getBytes());
        hapG.setScore(-3.0);
        final Haplotype hapGAnother = new Haplotype("GAAAG".getBytes());
        hapGAnother.setScore(-5.0);

        final Map<Allele, List<Haplotype>> alleleMapper = new LinkedHashMap<>();
        alleleMapper.put(ref, Arrays.asList(hapRef));
        alleleMapper.put(altC, Arrays.asList(hapC));
        alleleMapper.put(altT, Arrays.asList(hapT, hapTAnother));
        alleleMapper.put(altT2, Arrays.asList(hapT2, hapT2Another));
        alleleMapper.put(altG, Arrays.asList(hapG, hapGAnother));

        List<Allele> allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 5);
        Assert.assertEquals(allelesToKeep.size(), 5);

        Iterator<Allele> it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altC);
        Assert.assertEquals(it.next(), altT);
        Assert.assertEquals(it.next(), altT2);
        Assert.assertEquals(it.next(), altG);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 4);
        Assert.assertEquals(allelesToKeep.size(), 4);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altT);
        Assert.assertEquals(it.next(), altT2);
        Assert.assertEquals(it.next(), altG);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 3);
        Assert.assertEquals(allelesToKeep.size(), 3);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altT);
        Assert.assertEquals(it.next(), altT2);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 2);
        Assert.assertEquals(allelesToKeep.size(), 2);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altT);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 1);
        Assert.assertEquals(allelesToKeep.size(), 1);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
    }

    @Test
    public void testRemoveExcessiveAltAlleleFromVC(){
        final VariantContext originalVC = new VariantContextBuilder("source", "1", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("C", false), Allele.create("G", false))).make();

        final VariantContext reducedVC = HaplotypeCallerGenotypingEngine.removeExcessAltAllelesFromVC(originalVC, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("C", false)));

        Assert.assertEquals(reducedVC.getNAlleles(), 3);
        Assert.assertTrue(reducedVC.getAlleles().containsAll(Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("C", false))));
    }
}
