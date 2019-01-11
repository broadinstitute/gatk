package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.base.Strings;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.genotyper.HomogeneousPloidyModel;
import org.broadinstitute.hellbender.tools.walkers.genotyper.IndependentSampleGenotypesModel;
import org.broadinstitute.hellbender.tools.walkers.genotyper.PloidyModel;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;


public final class ReferenceConfidenceModelUnitTest extends GATKBaseTest {
    GenomeLocParser parser;
    final String RGID = "ID1";
    SAMReadGroupRecord rg;
    final String sample = "NA12878";
    final SampleList samples = SampleList.singletonSampleList(sample);
    SAMFileHeader header;
    ReferenceConfidenceModel model;

    @BeforeClass
    public void setUp() {
        header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);
        rg = new SAMReadGroupRecord(RGID);
        rg.setSample(sample);
        header.addReadGroup(rg);
        parser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @BeforeMethod
    public void setupModel() {
        model = new ReferenceConfidenceModel(samples, header, 10, -1);
    }

    @DataProvider(name = "CalcNIndelInformativeReadsData")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        { // very basic testing
            final String ref  = "ACGT";
            final String read = "ACGT";
            final String cigar = read.length() + "M";
            tests.add(new Object[]{read, cigar, ref, 1, Arrays.asList(1, 1, 1, 0)});
            tests.add(new Object[]{read, cigar, ref, 2, Arrays.asList(1, 1, 0, 0)});
            tests.add(new Object[]{read, cigar, ref, 3, Arrays.asList(1, 0, 0, 0)});
            tests.add(new Object[]{read, cigar, ref, 4, Arrays.asList(0, 0, 0, 0)});
        }

        { // actually interesting case where some sites aren't informative
            final String ref   = "TTAAAATT";
            final String read1 = "TTA";
            final String read2 = "TTAA";
            final String read3 = "TTAAA";
            final String read4 = "TTAAAA";
            final String read5 = "TTAAAAT";
            final String cigar1 = read1.length() + "M";
            final String cigar2 = read2.length() + "M";
            final String cigar3 = read3.length() + "M";
            final String cigar4 = read4.length() + "M";
            final String cigar5 = read5.length() + "M";

            tests.add(new Object[]{read1, cigar1, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read2, cigar2, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read3, cigar3, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read4, cigar4, ref, 1, Arrays.asList(1, 1, 0, 0, 0, 0, 0, 0)});
            tests.add(new Object[]{read5, cigar5, ref, 1, Arrays.asList(1, 1, 1, 1, 1, 1, 0, 0)});
        }

        for ( final String repeatUnit : Arrays.asList("A", "CA", "TAC", "TAGC", "TCAGA")) {
            final String anchor = Strings.repeat("G", repeatUnit.length());
            for ( int nUnits = 1; nUnits < 10; nUnits++ ) {
                final String repeat = Strings.repeat(repeatUnit, nUnits);
                final String ref = anchor + repeat + anchor;
                final String refCigar = ref.length() + "M";
                for ( int readLen = repeatUnit.length(); readLen < repeat.length(); readLen++ ) {
                    final String read = anchor + repeat.substring(0, readLen);
                    final String readCigar = read.length() + "M";
                    final List<Integer> expected = new LinkedList<>();
                    for ( int i = 0; i < anchor.length(); i++ ) expected.add(1);
                    for ( int i = 0; i < repeat.length(); i++ ) expected.add(readLen == repeat.length() ? 1 : 0);
                    for ( int i = 0; i < anchor.length(); i++ ) expected.add(0);
                    tests.add(new Object[]{read, readCigar, ref, repeatUnit.length(), expected});

                    final List<Integer> result = new ArrayList<>(Collections.nCopies(ref.length() - anchor.length(), 1));
                    result.addAll(Collections.nCopies(anchor.length(), 0));
                    tests.add(new Object[]{ref, refCigar, ref, repeatUnit.length(), result});
                }
            }
        }

        {//regression tests for an issue from late Aug 2018
            final String ref = "AATATCATCTTTGGTGTTT";
            final String read = "AATATCATTGGTG";
            final String cigar1 = read.length() + "M";
            final String cigar2 = "7M3D6M";
            //real issue is the informativeness for offset zero, but might as well run the rest of the offsets
            tests.add(new Object[]{read, cigar1, ref, 3, Arrays.asList(0,0,0,0,0,0,0,0,0,0,0,0,0)});
            tests.add(new Object[]{read, cigar2, ref, 3, Arrays.asList(1,1,1,1,1,1,0,0,0,0,0,0,0)});
        }
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "CalcNIndelInformativeReadsData")
    public void testCalcNIndelInformativeReads(final String readBases, final String cigar, final String ref, final int maxIndelSize, final List<Integer> expected ) {
        final byte qual = (byte)30;
        final byte[] quals = Utils.dupBytes(qual, readBases.length());

        for ( int i = 0; i < readBases.getBytes().length; i++ ) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(readBases.getBytes(), quals, cigar);
            final SimpleInterval loc = new SimpleInterval("20", i + 1, i + 1);
            final ReadPileup pileup = new ReadPileup(loc, Collections.singletonList(read), i);
            final int actual = model.calcNIndelInformativeReads(pileup, i, ref.getBytes(), maxIndelSize);
            Assert.assertEquals(actual, (int)expected.get(i), "failed at position " + i);
        }
    }

    @Test
    public void testWorstGL() {
        final GenotypeLikelihoods gq10 = GenotypeLikelihoods.fromPLField("0,10,100");
        final GenotypeLikelihoods gq20 = GenotypeLikelihoods.fromPLField("0,20,200");
        final GenotypeLikelihoods gq0 = GenotypeLikelihoods.fromPLField("20,0,200");

        Assert.assertSame(model.getGLwithWorstGQ(gq10, gq20), gq10);
        Assert.assertSame(model.getGLwithWorstGQ(gq20, gq10), gq10);
        Assert.assertSame(model.getGLwithWorstGQ(gq10, gq0), gq0);
        Assert.assertSame(model.getGLwithWorstGQ(gq0, gq10), gq0);
    }

    @Test
    public void testGetHeaderLines() throws Exception {
        final Set<VCFHeaderLine> vcfHeaderLines = model.getVCFHeaderLines();
        Assert.assertEquals(vcfHeaderLines.size(), 1);
        Assert.assertEquals(vcfHeaderLines.iterator().next(), new VCFSimpleHeaderLine(GATKVCFConstants.SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME, "Represents any possible alternative allele at this location"));
    }

    @Test
    public void testIndelLikelihoods() {
        GenotypeLikelihoods prev = model.getIndelPLs(HomoSapiensConstants.DEFAULT_PLOIDY,0);
        Assert.assertEquals(prev.getAsPLs(), new int[]{0, 0, 0});
        Assert.assertEquals(-10 * GenotypeLikelihoods.getGQLog10FromLikelihoods(GenotypeType.HOM_REF.ordinal()-1, prev.getAsVector()), 0.0);

        for ( int i = 1; i <= ReferenceConfidenceModel.MAX_N_INDEL_INFORMATIVE_READS; i++ ) {
            final GenotypeLikelihoods current = model.getIndelPLs(HomoSapiensConstants.DEFAULT_PLOIDY,i);
            final double prevGQ = -10 * GenotypeLikelihoods.getGQLog10FromLikelihoods(GenotypeType.HOM_REF.ordinal()-1, prev.getAsVector());
            final double currGQ = -10 * GenotypeLikelihoods.getGQLog10FromLikelihoods(GenotypeType.HOM_REF.ordinal()-1, current.getAsVector());
            Assert.assertTrue(prevGQ < currGQ, "GQ Failed with prev " + prev + " curr " + current + " at " + i);
            Assert.assertTrue(prev.getAsPLs()[1] < current.getAsPLs()[1], "het PL failed with prev " + prev + " curr " + current + " at " + i);
            Assert.assertTrue(prev.getAsPLs()[2] < current.getAsPLs()[2], "hom-var PL Failed with prev " + prev + " curr " + current + " at " + i);
//            logger.warn("result at " + i + " is " + current);
            prev = current;
        }
    }

    @Test
    public void testOverlappingVariantContext() {
        final VariantContext vc10 = GATKVariantContextUtils.makeFromAlleles("test", "1", 10, Arrays.asList("A", "C"));
        final VariantContext vc13 = GATKVariantContextUtils.makeFromAlleles("test", "1", 13, Arrays.asList("A", "C"));
        final VariantContext vc12_15 = GATKVariantContextUtils.makeFromAlleles("test", "1", 12, Arrays.asList("ACAT", "A"));
        final VariantContext vc18 = GATKVariantContextUtils.makeFromAlleles("test", "1", 18, Arrays.asList("A", "ACAT"));

        final List<VariantContext> calls = Arrays.asList(vc13, vc12_15, vc18, vc10);

        checkOverlapping(8, calls, null);
        checkOverlapping(9, calls, null);
        checkOverlapping(10, calls, vc10);
        checkOverlapping(11, calls, null);
        checkOverlapping(12, calls, vc12_15);
        checkOverlapping(13, calls, vc13);
        checkOverlapping(14, calls, vc12_15);
        checkOverlapping(15, calls, vc12_15);
        checkOverlapping(16, calls, null);
        checkOverlapping(17, calls, null);
        checkOverlapping(18, calls, vc18);
        checkOverlapping(19, calls, null);
        checkOverlapping(20, calls, null);
    }

    @DataProvider(name = "OffsetTestData")
    public Object[][] getOffsetTestData() {
        List<Object[]> tests = new ArrayList<>();

        final String read1 = "AATATCATTGGTG";
        final String cigar1 = "13M";
        final String cigar2 = "7M3D6M";
        final String cigar3 = "7M3I3M";
        final String cigar4 = "7M3D3I3M";
        //offset inside soft clip is not allowed
        tests.add(new Object[]{read1, cigar1, 5, 5});
        tests.add(new Object[]{read1, cigar2, 5, 5});
        tests.add(new Object[]{read1, cigar2, 11, 14});
        tests.add(new Object[]{read1, cigar3, 5, 5});
        tests.add(new Object[]{read1, cigar3, 11, 8});
        tests.add(new Object[]{read1, cigar4, 5, 5});
        tests.add(new Object[]{read1, cigar4, 11, 11});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "OffsetTestData")
    public void testCigarModifiedOffset(final String readBases, final String cigar, final int pileupOffset, final int expectedNewOffset) {
        final byte qual = (byte)30;
        final byte[] quals = Utils.dupBytes(qual, readBases.length());

        final GATKRead read = ArtificialReadUtils.createArtificialRead(readBases.getBytes(), quals, cigar);
        final PileupElement pe = PileupElement.createPileupForReadAndOffset(read, pileupOffset);
        final int newOffset = model.getCigarModifiedOffset(pe);
        Assert.assertEquals(newOffset, expectedNewOffset);
    }

    private void checkOverlapping(final int pos, Collection<VariantContext> calls, final VariantContext expected) {
        final GenomeLoc loc = parser.createGenomeLoc(parser.getSequenceDictionary().getSequences().get(0).getSequenceName(), pos, pos);
        final VariantContext actual = GATKVariantContextUtils.getOverlappingVariantContext(loc, calls);
        Assert.assertEquals(actual, expected);
    }

    //
    // test reference calculation
    //
    private class RefConfData {
        final String ref;
        final int extension;
        final Haplotype refHap;
        final SimpleInterval refLoc, paddedRefLoc;
        final AssemblyRegion region;
        int readCounter = 0;

        private RefConfData(String ref, int extension) {
            this.ref = ref;
            this.extension = extension;

            refLoc = new SimpleInterval("1", getStart(), getEnd());
            paddedRefLoc = new SimpleInterval("1", getStart() - extension, getEnd() + extension);
            region = new AssemblyRegion(getRefLoc(), extension, header);
            final String pad = Strings.repeat("N", extension);
            refHap = ReferenceConfidenceModel.createReferenceHaplotype(getActiveRegion(), (pad + ref + pad).getBytes(), getPaddedRefLoc());
        }

        public SimpleInterval getRefLoc() { return refLoc; }
        public SimpleInterval getPaddedRefLoc() { return paddedRefLoc; }
        public AssemblyRegion getActiveRegion() { return region; }
        public Haplotype getRefHap() { return refHap; }
        public int getStart() { return 100; }
        public int getEnd() { return getStart() + getRefLength() - 1; }
        public byte[] getRefBases() { return ref.getBytes(); }
        public int getRefLength() { return ref.length(); }

        public GATKRead makeRead(final int start, final int length) {
            final byte[] quals = Utils.dupBytes((byte)30, length);
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read " + readCounter++, 0, start + getStart(), ref.substring(start, start + length).getBytes(), quals, length + "M");
            read.setReadGroup(rg.getId());
            return read;
        }
    }


    @DataProvider(name = "RefConfidenceData")
    public Object[][] makeRefConfidenceData() {
        List<Object[]> tests = new ArrayList<>();

        for ( int i = 0; i < 10; i++ ) {
            for ( final int extension : Arrays.asList(0, 10) ) {
                tests.add(new Object[]{i, extension});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RefConfidenceData")
    public void testRefConfidenceBasic(final int nReads, final int extension) {
        final RefConfData data = new RefConfData("ACGTAACCGGTT", extension);
        final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
        final List<VariantContext> calls = Collections.emptyList();

        for ( int i = 0; i < nReads; i++ ) {
            data.getActiveRegion().add(data.makeRead(0, data.getRefLength()));
        }

        final ReadLikelihoods<Haplotype> likelihoods = createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final IndependentSampleGenotypesModel genotypingModel = new IndependentSampleGenotypesModel();
        final List<Integer> expectedDPs = Collections.nCopies(data.getActiveRegion().getSpan().size(), nReads);
        final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, calls, false, Collections.emptyList());
        checkReferenceModelResult(data, contexts, expectedDPs, calls);
    }

    @Test
    public void testRefConfidencePartialReads() {

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final IndependentSampleGenotypesModel genotypingModel = new IndependentSampleGenotypesModel();
        final String ref = "ACGTAACCGGTT";
        for ( int readLen = 3; readLen < ref.length(); readLen++ ) {
            for ( int start = 0; start < ref.length() - readLen; start++ ) {
                final RefConfData data = new RefConfData(ref, 0);
                final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
                final List<VariantContext> calls = Collections.emptyList();

                data.getActiveRegion().add(data.makeRead(start, readLen));
                final ReadLikelihoods<Haplotype> likelihoods = createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

                final List<Integer> expectedDPs = new ArrayList<>(Collections.nCopies(data.getActiveRegion().getSpan().size(), 0));
                for ( int i = start; i < readLen + start; i++ ) expectedDPs.set(i, 1);
                final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, calls);
                checkReferenceModelResult(data, contexts, expectedDPs, calls);
            }
        }
    }

    @Test
    public void testRefConfidenceWithCalls() {
        final RefConfData xxxdata = new RefConfData("ACGTAACCGGTT", 0);
        final int start = xxxdata.getStart();
        final int stop = xxxdata.getEnd();

        final PloidyModel ploidyModel = new HomogeneousPloidyModel(samples,2);
        final IndependentSampleGenotypesModel genotypingModel = new IndependentSampleGenotypesModel();

        for ( int nReads = 0; nReads < 2; nReads++ ) {

            final VariantContext vcStart = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start, Arrays.asList("A", "C"));
            final VariantContext vcEnd = GATKVariantContextUtils.makeFromAlleles("test", "chr1", stop, Arrays.asList("A", "C"));
            final VariantContext vcMiddle = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 2, Arrays.asList("A", "C"));
            final VariantContext vcDel = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 4, Arrays.asList("AAC", "A"));
            final VariantContext vcIns = GATKVariantContextUtils.makeFromAlleles("test", "chr1", start + 8, Arrays.asList("G", "GCG"));

            final List<VariantContext> allCalls = Arrays.asList(vcStart, vcEnd, vcMiddle, vcDel, vcIns);

            for ( int n = 1; n <= allCalls.size(); n++ ) {
                for ( final List<VariantContext> calls : Utils.makePermutations(allCalls, n, false) ) {
                    final RefConfData data = new RefConfData("ACGTAACCGGTT", 0);
                    final List<Haplotype> haplotypes = Arrays.asList(data.getRefHap());
                    for ( int i = 0; i < nReads; i++ ) {
                        data.getActiveRegion().add(data.makeRead(0, data.getRefLength()));
                    }

                    final ReadLikelihoods<Haplotype> likelihoods = createDummyStratifiedReadMap(data.getRefHap(), samples, data.getActiveRegion());

                    final List<Integer> expectedDPs = Collections.nCopies(data.getActiveRegion().getSpan().size(), nReads);
                    final List<VariantContext> contexts = model.calculateRefConfidence(data.getRefHap(), haplotypes, data.getPaddedRefLoc(), data.getActiveRegion(), likelihoods, ploidyModel, calls);
                    checkReferenceModelResult(data, contexts, expectedDPs, calls);
                }
            }
        }
    }

    /**
     * Create a context that maps each read to the reference haplotype with log10 L of 0
     * @param refHaplotype a non-null reference haplotype
     * @param samples a list of all samples
     * @param region the active region containing reads
     * @return a map from sample -> PerReadAlleleLikelihoodMap that maps each read to ref
     */
    public static ReadLikelihoods<Haplotype> createDummyStratifiedReadMap(final Haplotype refHaplotype,
                                                                          final SampleList samples,
                                                                          final AssemblyRegion region) {
        return new ReadLikelihoods<>(samples, new IndexedAlleleList<>(refHaplotype),
                splitReadsBySample(samples, region.getReads(), region.getHeader()));
    }

    public static Map<String, List<GATKRead>> splitReadsBySample( final SampleList samplesList, final Collection<GATKRead> reads , final SAMFileHeader header) {
        final Map<String, List<GATKRead>> returnMap = new HashMap<>();
        final int sampleCount = samplesList.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            returnMap.put(samplesList.getSample(i), new ArrayList<>());
        }

        for( final GATKRead read : reads ) {
            returnMap.get(ReadUtils.getSampleName(read, header)).add(read);
        }

        return returnMap;
    }

    private void checkReferenceModelResult(final RefConfData data, final List<VariantContext> contexts, final List<Integer> expectedDPs, final List<VariantContext> calls) {
        Assert.assertNotNull(contexts);

        final SimpleInterval loc = data.getActiveRegion().getExtendedSpan();
        final List<Boolean> seenBP = new ArrayList<>(Collections.nCopies(data.getActiveRegion().getSpan().size(), false));

        for ( int i = 0; i < loc.size(); i++ ) {
            final GenomeLoc curPos = parser.createGenomeLoc(loc.getContig(), loc.getStart() + i);
            final VariantContext call = GATKVariantContextUtils.getOverlappingVariantContext(curPos, calls);
            final VariantContext refModel = GATKVariantContextUtils.getOverlappingVariantContext(curPos, contexts);

            if ( ! data.getActiveRegion().getSpan().contains(curPos) ) {
                // part of the extended interval, but not the full interval
                Assert.assertNull(refModel);
                continue;
            }

            if ( call != null ) {
                if (call.isVariant() && refModel.getType() ==  VariantContext.Type.SYMBOLIC ) {
                    //Assert.assertEquals(refModel, call, "Should have found call " + call + " but found " + refModel + " instead");
                    Assert.assertTrue(call.getReference().length() > 1); // must be a deletion.
                    Assert.assertTrue(call.getStart() < refModel.getStart()); // the deletion must not start at the same position
                    Assert.assertEquals(call.getReference().getBaseString().substring(refModel.getStart() - call.getStart(),
                                refModel.getStart() - call.getStart() + 1), refModel.getReference().getBaseString(), "" + data.getRefHap()); // the reference must be the same.
                    Assert.assertTrue(refModel.getGenotype(0).getGQ() <= 0); // No confidence in the reference hom-ref call across the deletion
                    Assert.assertEquals(refModel.getAlleles().size(),2); // the reference and the lonelly <NON_REF>
                    Assert.assertEquals(refModel.getAlleles().get(1), Allele.NON_REF_ALLELE);
                } else {
                    Assert.assertEquals(refModel, call, "Should have found call " + call + " but found " + refModel + " instead");
                }

            } else {
                final int expectedDP = expectedDPs.get(curPos.getStart() - data.getActiveRegion().getSpan().getStart());
                Assert.assertEquals(refModel.getStart(), loc.getStart() + i);
                Assert.assertEquals(refModel.getEnd(), loc.getStart() + i);
                Assert.assertFalse(refModel.hasLog10PError());
                Assert.assertEquals(refModel.getAlternateAlleles().size(), 1);
                Assert.assertEquals(refModel.getAlternateAllele(0), Allele.NON_REF_ALLELE);
                Assert.assertTrue(refModel.hasGenotype(sample));

                final Genotype g = refModel.getGenotype(sample);
                Assert.assertTrue(g.hasAD());
                Assert.assertTrue(g.hasDP());
                Assert.assertEquals(g.getDP(), expectedDP);
                Assert.assertTrue(g.hasGQ());
                Assert.assertTrue(g.hasPL());
            }

            final VariantContext vc = call == null ? refModel : call;
            if ( curPos.getStart() == vc.getStart() ) {
                for ( int pos = vc.getStart(); pos <= vc.getEnd(); pos++ ) {
                    final int j = pos - data.getActiveRegion().getSpan().getStart();
                    Assert.assertFalse(seenBP.get(j));
                    seenBP.set(j, true);
                }
            }
        }

        for ( int i = 0; i < seenBP.size(); i++ ) {
            Assert.assertEquals((boolean)seenBP.get(i), true);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRefVsAnyResultNotNegative() throws Exception {
        new RefVsAnyResult(-1);
    }
    @Test
    public void testRefVsAnyResultConstructor() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        Assert.assertEquals(res.getAD().length, 2);
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood().length, 3);
        Assert.assertEquals(res.getDP(), 0);
        Assert.assertEquals(res.getAD(), new int[]{0, 0});
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{0.0, 0.0, 0.0});
    }

    @Test
    public void testRefVsAnyResultADInc() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.refDepth += 2;
        res.nonRefDepth += 3;
        Assert.assertEquals(res.getAD(), new int[]{2, 3});
    }

    @Test
    public void testRefVsAnyResultCapByHomRefLikelihood() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.genotypeLikelihoods[0] += 100;
        res.genotypeLikelihoods[1] += 200;
        res.genotypeLikelihoods[2] += 60;
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{100.0, 100.0, 60.0});
    }

    @Test
    public void testRefVsAnyResultArrays() throws Exception {
        final RefVsAnyResult res = new RefVsAnyResult(3);
        res.refDepth += 2;
        res.nonRefDepth += 3;
        final int[] adArray = res.getAD();
        Assert.assertEquals(adArray, new int[]{2, 3});

        adArray[0] = 17;
        Assert.assertEquals(res.getAD(), new int[]{2, 3}); //verify that the ad array is a copy

        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{0, 0, 0});
        double[] liks = res.getGenotypeLikelihoodsCappedByHomRefLikelihood();
        liks[0] = 19;
        Assert.assertEquals(res.getGenotypeLikelihoodsCappedByHomRefLikelihood(), new double[]{0, 0, 0}); //verify that the GL array is a copy
    }

}
