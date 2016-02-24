package org.broadinstitute.hellbender.utils.genotyper;

import com.google.common.base.Strings;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public final class PerReadAlleleLikelihoodMapUnitTest extends BaseTest {

    // example fasta index file, can be deleted if you don't use the reference
    private IndexedFastaSequenceFile seq;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg19MiniReference));
    }

    @Test()
    public void testMultiAlleleWithHomLiks_2() {
        final ReadPileup pileup = makeArtificialPileup("first");

        Allele base_A = Allele.create(BaseUtils.Base.A.base);
        Allele base_C = Allele.create(BaseUtils.Base.C.base);
        Allele base_T = Allele.create(BaseUtils.Base.T.base);

        final List<Allele> allAlleles = Arrays.asList(base_A, base_C, base_T);
        final double likA = -0.04;
        final double likNotA = -3.0;
        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        for ( final PileupElement e : pileup ) {
            for ( final Allele allele : allAlleles ) {
                Double likelihood = allele == base_A ? likA : likNotA;
                perReadAlleleLikelihoodMap.add(e, allele, likelihood);
            }
        }

        Assert.assertEquals(perReadAlleleLikelihoodMap.getAllelesSet(), new HashSet<>(allAlleles));
        Assert.assertEquals(perReadAlleleLikelihoodMap.size(), pileup.size());
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().keySet().size(), 3);
        Map<Allele,List<GATKRead>> shouldBeAllA = perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap();
        Assert.assertEquals(shouldBeAllA.get(base_A).size(), pileup.size());
        Assert.assertEquals(shouldBeAllA.get(base_C).size(), 0);
        Assert.assertEquals(shouldBeAllA.get(base_T).size(), 0);

        final Map<GATKRead, Map<Allele, Double>> readMap = perReadAlleleLikelihoodMap.getLikelihoodReadMap();
        Assert.assertEquals(readMap.size(), pileup.size());

        final Map<Allele, Double> likMap = perReadAlleleLikelihoodMap.getLikelihoods(pileup.iterator().next());
        Assert.assertEquals(likMap.size(), 3); //three allele
        Assert.assertEquals(likMap.get(base_A), likA);
        Assert.assertEquals(likMap.get(base_C), likNotA);

        ReadPileup newPileup = makeArtificialPileup("second");
        Assert.assertNull(perReadAlleleLikelihoodMap.getLikelihoods(newPileup.iterator().next()));

        Assert.assertNotNull(perReadAlleleLikelihoodMap.toString()); //checking blowup

        perReadAlleleLikelihoodMap.clear();
        Assert.assertEquals(perReadAlleleLikelihoodMap.size(), 0);

        Assert.assertNotNull(perReadAlleleLikelihoodMap.toString()); //checking blowup
    }

    @Test()
    public void testMultiAlleleWithHomLiks() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final Locatable myLocation = new SimpleInterval("1", 10, 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKRead> reads = new LinkedList<>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte) 'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other read the first base is a C

            // set the read's bases and quals
            read.setBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        final ReadPileup pileup = new ReadPileup(myLocation, reads, 0);
        Allele base_A = Allele.create(BaseUtils.Base.A.base);
        Allele base_C = Allele.create(BaseUtils.Base.C.base);
        Allele base_T = Allele.create(BaseUtils.Base.T.base);

        final List<Allele> allAlleles = Arrays.asList(base_A, base_C, base_T);
        final double likA = -0.04;
        final double likNotA = -3.0;
        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        for ( final PileupElement e : pileup ) {
            for ( final Allele allele : allAlleles ) {
                Double likelihood = allele == base_A ? likA : likNotA;
                perReadAlleleLikelihoodMap.add(e, allele, likelihood);
            }
        }

        Assert.assertEquals(perReadAlleleLikelihoodMap.getAllelesSet(), new HashSet<>(allAlleles));
        Assert.assertEquals(perReadAlleleLikelihoodMap.size(), pileup.size());
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().keySet().size(), 3);
        Map<Allele,List<GATKRead>> shouldBeAllA = perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap();
        Assert.assertEquals(shouldBeAllA.get(base_A).size(), pileup.size());
        Assert.assertEquals(shouldBeAllA.get(base_C).size(), 0);
        Assert.assertEquals(shouldBeAllA.get(base_T).size(), 0);

        //getLikelihoodReadMap
        final Map<GATKRead, Map<Allele, Double>> likelihoodReadMap = perReadAlleleLikelihoodMap.getLikelihoodReadMap();
        Assert.assertEquals(likelihoodReadMap.size(), pileupSize);

        Assert.assertNotNull(perReadAlleleLikelihoodMap.toString());//just checking non-null results. The contents are only human readible.

        final Set<Allele> allelesSet = perReadAlleleLikelihoodMap.getAllelesSet();
        Assert.assertEquals(allelesSet, new HashSet<>(Arrays.asList(base_A, base_C, base_T)));

        final int oldSize = perReadAlleleLikelihoodMap.size();
        perReadAlleleLikelihoodMap.performPerAlleleDownsampling(-0.5); //no-op
        final int newSize = perReadAlleleLikelihoodMap.size();
        Assert.assertEquals(oldSize, newSize);

        final MostLikelyAllele mla = PerReadAlleleLikelihoodMap.getMostLikelyAllele(Collections.emptyMap());
        Assert.assertTrue(mla.getMostLikelyAllele().isNoCall());
        Assert.assertNull(mla.getSecondMostLikelyAllele());
        Assert.assertEquals(mla.getLog10LikelihoodOfMostLikely(), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(mla.getLog10LikelihoodOfSecondBest(), Double.NEGATIVE_INFINITY);

        final MostLikelyAllele mla2 = PerReadAlleleLikelihoodMap.getMostLikelyAllele(Collections.singletonMap(base_A, -0.04));
        Assert.assertFalse(mla2.getMostLikelyAllele().isNoCall());
        Assert.assertNull(mla2.getSecondMostLikelyAllele());
        Assert.assertEquals(mla2.getLog10LikelihoodOfMostLikely(), -0.04);
        Assert.assertEquals(mla2.getLog10LikelihoodOfSecondBest(), Double.NEGATIVE_INFINITY);

        //clear makes it empty
        perReadAlleleLikelihoodMap.clear();
        Assert.assertTrue(perReadAlleleLikelihoodMap.isEmpty());
        Assert.assertNull(perReadAlleleLikelihoodMap.getMostLikelyDiploidAlleles());
        Assert.assertEquals(perReadAlleleLikelihoodMap.size(), 0);
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().keySet().size(), 0);

        final PileupElement pu = pileup.iterator().next();
        final Map<Allele, Double> ll1 = perReadAlleleLikelihoodMap.getLikelihoodsAssociatedWithPileupElement(pu);
        Assert.assertNull(ll1);

        final Map<Allele, Double> ll2 = perReadAlleleLikelihoodMap.getLikelihoods(pu);
        Assert.assertNull(ll2);

        Assert.assertNotNull(perReadAlleleLikelihoodMap.toString());//just checking non-null results. The contents are only human readible.
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testNoGoodLikelihood() throws Exception {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final Locatable myLocation = new SimpleInterval("1", 10, 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKRead> reads = new LinkedList<>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte) 'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other read the first base is a C

            // set the read's bases and quals
            read.setBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        final ReadPileup pileup = new ReadPileup(myLocation, reads, 0);
        Allele base_A = Allele.create(BaseUtils.Base.A.base);
        Allele base_C = Allele.create(BaseUtils.Base.C.base);
        Allele base_T = Allele.create(BaseUtils.Base.T.base);

        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        for ( final PileupElement e : pileup ) {
            for ( final Allele allele : Arrays.asList(base_A, base_C, base_T) ) {
                double likelihood = Double.NEGATIVE_INFINITY;
                perReadAlleleLikelihoodMap.add(e,allele,likelihood);
            }
        }
        perReadAlleleLikelihoodMap.getMostLikelyDiploidAlleles(); //this should blow up because all likelihoods are infinities
    }

    private ReadPileup makeArtificialPileup(final String prefix) {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final Locatable myLocation = new SimpleInterval("1", 10, 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKRead> reads = new LinkedList<>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, prefix + "_myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte) 'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other read the first base is a C

            // set the read's bases and quals
            read.setBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        return new ReadPileup(myLocation, reads, 0);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadLikelihood(){
        final ReadPileup pileup = makeArtificialPileup("first");
        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        Allele base_A = Allele.create(BaseUtils.Base.A.base);
        final Allele allele = base_A;
        PileupElement e = pileup.iterator().next();
        perReadAlleleLikelihoodMap.add(e,allele, 0.1);  //boom: positive likelihood
    }

    @Test()
    public void testMultiAlleleWithHetLiks() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final Locatable myLocation = new SimpleInterval("1", 10, 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKRead> reads = new LinkedList<>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte)'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other base is a C

            // set the read's bases and quals
            read.setBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        final ReadPileup pileup = new ReadPileup(myLocation, reads, 0);
        Allele base_A = Allele.create(BaseUtils.Base.A.base);
        Allele base_C = Allele.create(BaseUtils.Base.C.base);
        Allele base_T = Allele.create(BaseUtils.Base.T.base);

        PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = new PerReadAlleleLikelihoodMap();
        int idx = 0;
        for ( final PileupElement e : pileup ) {
            for ( final Allele allele : Arrays.asList(base_A, base_C, base_T) ) {
                Double likelihood;
                if ( idx % 2 == 0 )
                    likelihood = allele == base_A ? -0.04 : -3.0;
                else
                    likelihood = allele == base_C ? -0.04 : -3.0;
                perReadAlleleLikelihoodMap.add(e,allele,likelihood);
            }
            idx++;
        }

        Assert.assertEquals(perReadAlleleLikelihoodMap.size(), pileup.size());
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().keySet().size(), 3);
        Map<Allele,List<GATKRead>> halfAhalfC = perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap();
        Assert.assertEquals(halfAhalfC.get(base_A).size(), pileup.size() / 2);
        Assert.assertEquals(halfAhalfC.get(base_C).size(), pileup.size() / 2);
        Assert.assertEquals(halfAhalfC.get(base_T).size(), 0);

        // make sure the likelihoods are retrievable

        idx = 0;
        for ( final PileupElement e : pileup ) {
            Assert.assertTrue(perReadAlleleLikelihoodMap.containsPileupElement(e));
            Map<Allele,Double> likelihoods = perReadAlleleLikelihoodMap.getLikelihoodsAssociatedWithPileupElement(e);
            for ( final Allele allele : Arrays.asList(base_A, base_C, base_T) ) {
                Double expLik;
                if ( idx % 2 == 0 )
                    expLik = allele == base_A ? -0.04 : -3.0;
                else
                    expLik = allele == base_C ? -0.04 : -3.0;
                Assert.assertEquals(likelihoods.get(allele), expLik);
            }
            idx++;
        }

        // and test downsampling for good measure

        final List<GATKRead> excessReads = new LinkedList<>();
        int prevSize = perReadAlleleLikelihoodMap.size();
        for ( int i = 0; i < 10 ; i++ ) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myExcessRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte)'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other base is a C

            // set the read's bases and quals
            read.setBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            for ( final Allele allele : Arrays.asList(base_A, base_C, base_T) ) {
                perReadAlleleLikelihoodMap.add(read,allele,allele==base_A ? -0.04 : -3.0);
            }
            Assert.assertEquals(perReadAlleleLikelihoodMap.size(), 1 + prevSize);
            prevSize = perReadAlleleLikelihoodMap.size();
        }

        Assert.assertEquals(perReadAlleleLikelihoodMap.size(), pileup.size() + 10);
        Assert.assertEquals(perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap().get(base_A).size(), 60);
        perReadAlleleLikelihoodMap.performPerAlleleDownsampling(0.1);
        Assert.assertEquals(perReadAlleleLikelihoodMap.size(), (int) (0.9 * (pileup.size() + 10)));

        Map<Allele,List<GATKRead>> downsampledStrat = perReadAlleleLikelihoodMap.getAlleleStratifiedReadMap();
        Assert.assertEquals(downsampledStrat.get(base_A).size(), (pileup.size() / 2) - 1);
        Assert.assertEquals(downsampledStrat.get(base_C).size(), (pileup.size() / 2));
        Assert.assertEquals(downsampledStrat.get(base_T).size(), 0);

        perReadAlleleLikelihoodMap.performPerAlleleDownsampling(1.2);   //effectively clear
        Assert.assertTrue(perReadAlleleLikelihoodMap.isEmpty());
    }

    @DataProvider(name = "PoorlyModelledReadData")
    public Object[][] makePoorlyModelledReadData() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        tests.add(new Object[]{10, 0.1, false, Arrays.asList(0.0)});
        tests.add(new Object[]{10, 0.1, true, Arrays.asList(-10.0)});
        tests.add(new Object[]{10, 0.1, false, Arrays.asList(0.0, -10.0)});
        tests.add(new Object[]{10, 0.1, true, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.1, false, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.01, true, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.01, false, Arrays.asList(-5.0, -10.0, -3.0)});
        tests.add(new Object[]{100, 0.01, false, Arrays.asList(-5.0, -10.0, -2.0)});
        tests.add(new Object[]{100, 0.01, true, Arrays.asList(-5.0, -10.0, -4.2)});
        tests.add(new Object[]{100, 0.001, true, Arrays.asList(-5.0, -10.0)});
        tests.add(new Object[]{100, 0.001, false, Arrays.asList(-5.0, -10.0, 0.0)});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "PoorlyModelledReadData")
    public void testPoorlyModelledRead(final int readLen, final double maxErrorRatePerBase, final boolean expected, final List<Double> log10likelihoods) {
        final byte[] bases = Utils.dupBytes((byte)'A', readLen);
        final byte[] quals = Utils.dupBytes((byte) 40, readLen);

        final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, readLen + "M");

        final boolean actual = PerReadAlleleLikelihoodMap.readIsPoorlyModelled(read, log10likelihoods, maxErrorRatePerBase);
        Assert.assertEquals(actual, expected);
    }


    @DataProvider(name = "RemovingPoorlyModelledReadData")
    public Object[][] makeRemovingPoorlyModelledReadData() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        final int readLen = 10;
        for ( int nReads = 0; nReads < 4; nReads++ ) {
            for ( int nBad = 0; nBad <= nReads; nBad++ ) {
                final int nGood = nReads - nBad;
                tests.add(new Object[]{readLen, nReads, nBad, nGood});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RemovingPoorlyModelledReadData")
    public void testRemovingPoorlyModelledReads(final int readLen, final int nReads, final int nBad, final int nGood) {
        final PerReadAlleleLikelihoodMap map = new PerReadAlleleLikelihoodMap();
        final Set<GATKRead> goodReads = new HashSet<>();
        final Set<GATKRead> badReads = new HashSet<>();
        for ( int readI = 0; readI < nReads; readI++ ) {
            final boolean bad = readI < nBad;
            final double likelihood = bad ? -100.0 : 0.0;

            final byte[] bases = Utils.dupBytes((byte)'A', readLen);
            final byte[] quals = Utils.dupBytes((byte) 40, readLen);

            final Allele allele = Allele.create(Strings.repeat("A", readI + 1));

            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, readLen + "M");
            read.setName("readName" + readI);
            map.add(read, allele, likelihood);
            (bad ? badReads : goodReads).add(read);
        }

        final List<GATKRead> removedReads = map.filterPoorlyModelledReads(0.01);
        Assert.assertEquals(removedReads.size(), nBad, "nBad " + nBad + " nGood " + nGood);
        Assert.assertEquals(new HashSet<>(removedReads), badReads, "nBad " + nBad + " nGood " + nGood);
        Assert.assertEquals(map.size(), nGood, "nBad " + nBad + " nGood " + nGood);
        Assert.assertTrue(map.getReads().containsAll(goodReads), "nBad " + nBad + " nGood " + nGood);
        Assert.assertEquals(map.getReads().size(), nGood, "nBad " + nBad + " nGood " + nGood);
    }

    @DataProvider(name = "MostLikelyAlleleData")
    public Object[][] makeMostLikelyAlleleData() {
        List<Object[]> tests = new ArrayList<>();

        final Allele a = Allele.create("A");
        final Allele c = Allele.create("C");
        final Allele g = Allele.create("G");

        tests.add(new Object[]{Arrays.asList(a), Arrays.asList(Arrays.asList(0.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c), Arrays.asList(Arrays.asList(0.0, -1.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c), Arrays.asList(Arrays.asList(-1.0, 0.0)), c, c});
        tests.add(new Object[]{Arrays.asList(a, c, g), Arrays.asList(Arrays.asList(0.0, 0.0, -10.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c, g), Arrays.asList(Arrays.asList(0.0, 0.0, -10.0)), a, a});
        tests.add(new Object[]{Arrays.asList(a, c, g),
                Arrays.asList(
                        Arrays.asList(0.0, -10.0, -10.0),
                        Arrays.asList(-100.0, 0.0, -10.0)),
                c, a});
        tests.add(new Object[]{Arrays.asList(a, c, g),
                Arrays.asList(
                        Arrays.asList(0.0, -10.0, -10.0),
                        Arrays.asList(-20.0, 0.0, -100.0)),
                c, a});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MostLikelyAlleleData")
    public void testMostLikelyAllele(final List<Allele> alleles, final List<List<Double>> perReadlikelihoods, final Allele best, final Allele second) {
        final PerReadAlleleLikelihoodMap map = new PerReadAlleleLikelihoodMap();

        for ( int readI = 0; readI < perReadlikelihoods.size(); readI++ ) {
            final List<Double> likelihoods = perReadlikelihoods.get(readI);

            final byte[] bases = Utils.dupBytes((byte)'A', 10);
            final byte[] quals = Utils.dupBytes((byte) 30, 10);
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, "10M");
            read.setName("readName" + readI);

            for ( int i = 0; i < alleles.size(); i++ ) {
                final Allele allele = alleles.get(i);
                final double likelihood = likelihoods.get(i);
                map.add(read, allele, likelihood);
            }
        }

        final Pair<Allele, Allele> mla = getMostLikelyDiploidAlleles(map);
        Assert.assertEquals(mla.getLeft(), best);
        Assert.assertEquals(mla.getRight(), second);
    }

    /**
     * Get the most likely alleles estimated across all reads in this object
     *
     * Takes the most likely two alleles according to their diploid genotype likelihoods.  That is, for
     * each allele i and j we compute p(D | i,j) where D is the read likelihoods.  We track the maximum
     * i,j likelihood and return an object that contains the alleles i and j as well as the max likelihood.
     *
     * Note that the second most likely diploid genotype is not tracked so the resulting MostLikelyAllele
     * doesn't have a meaningful get best likelihood.
     *
     * @return a pair of Alleles, first is most likely second is second most likely. Can be the same ones.
     */
    private static Pair<Allele, Allele> getMostLikelyDiploidAlleles(final PerReadAlleleLikelihoodMap map) {
        if ( map.isEmpty() ) {
            return null;
        }

        Allele bestAllele1= null;
        Allele bestAllele2= null;
        double maxElement = Double.NEGATIVE_INFINITY;
        for( int i = 0; i < map.alleleCount(); i++ ) {
            final Allele allele_i = map.getAllele(i);
            for( int j = 0; j <= i; j++ ) {
                final Allele allele_j = map.getAllele(j);

                double haplotypeLikelihood = 0.0;
                for( final Map<Allele,Double> alleleMap : map.getAlleleToLikelihoodMaps() ) {
                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    final double log10_likelihood_i = alleleMap.get(allele_i);
                    final double log10_likelihood_j = alleleMap.get(allele_j);
                    haplotypeLikelihood += MathUtils.approximateLog10SumLog10(log10_likelihood_i, log10_likelihood_j) + MathUtils.LOG10_ONE_HALF;

                    // fast exit.  If this diploid pair is already worse than the max, just stop and look at the next pair
                    if ( haplotypeLikelihood < maxElement ) {
                        break;
                    }
                }

                // keep track of the max element and associated alleles
                if ( haplotypeLikelihood > maxElement ) {
                    bestAllele1 = allele_i;
                    bestAllele2 = allele_j;
                    maxElement = haplotypeLikelihood;
                }
            }
        }

        if ( maxElement == Double.NEGATIVE_INFINITY ) {
            throw new IllegalStateException("max likelihood is " + maxElement + " indicating something has gone wrong");
        }
        if (bestAllele1 == null || bestAllele2 == null){
            throw new IllegalStateException("best alleles are null, indicating something has gone wrong");
        }

        return Pair.of(bestAllele1, bestAllele2);
    }


    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIllegalLikelihood() throws Exception {
        final PerReadAlleleLikelihoodMap map = new PerReadAlleleLikelihoodMap();
        final byte[] bases = Utils.dupBytes((byte)'A', 10);
        final byte[] quals = Utils.dupBytes((byte) 30, 10);
        final int readLen = 100;
        final GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, readLen + "M");

        final Allele allele = Allele.create("A");
        map.add(read, allele, 3.0);
    }

}
