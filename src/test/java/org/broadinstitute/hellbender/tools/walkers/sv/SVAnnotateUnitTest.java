package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.tribble.bed.FullBEDFeature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.compress.utils.Sets;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.ClosedSVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfGeneFeature;
import org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfTranscriptFeature;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

import static java.util.Objects.isNull;

public class SVAnnotateUnitTest extends GATKBaseTest {
    @DataProvider(name="variantFeatureComparisons")
    public Object[][] getVariantFeatureComparisonsData() {
        return new Object[][] {
                { new SimpleInterval("chr1", 1, 2000), new SimpleInterval("chr2", 100, 200), false, 0 },
                { new SimpleInterval("chrX", 6000, 8000), new SimpleInterval("chrX", 30, 400), false, 0 },
                { new SimpleInterval("chr3", 10000, 30000), new SimpleInterval("chr3", 9000, 20000), false, 1 },
                { new SimpleInterval("chrY", 500, 700), new SimpleInterval("chrY", 600, 800), false, 1 },
                { new SimpleInterval("chr22", 40, 90), new SimpleInterval("chr22", 25, 700), false, 2 },
                { new SimpleInterval("chr19", 33, 33), new SimpleInterval("chr19", 22, 44), false, 2 },
                { new SimpleInterval("chr10", 900, 4000), new SimpleInterval("chr10", 2000, 2200), true, 0 },
        };
    }

    @Test(dataProvider = "variantFeatureComparisons")
    public void testVariantFeatureComparisons(
            final SimpleInterval variantInterval,
            final SimpleInterval featureInterval,
            final boolean expectedVariantSpansFeature,
            final int expectedNumBreakpointsInsideFeature)
    {
        final boolean actualVariantSpansFeature = SVAnnotate.variantSpansFeature(variantInterval, featureInterval);
        Assert.assertEquals(expectedVariantSpansFeature, actualVariantSpansFeature);

        final int actualNumBreakpointsInsideFeature = SVAnnotate.countBreakendsInsideFeature(variantInterval, featureInterval);
        Assert.assertEquals(expectedNumBreakpointsInsideFeature, actualNumBreakpointsInsideFeature);
    }

    // transcription start sites and initTree() method for testing annotateNearestTranscriptionStartSite()
    private static final ClosedSVInterval[] transcriptionStartSites = {
            new ClosedSVInterval(0, 100, 101),
            new ClosedSVInterval(1, 150, 151),
            new ClosedSVInterval(1, 200, 201),
            new ClosedSVInterval(1, 250, 251),
            new ClosedSVInterval(2, 1, 2)
    };

    private static SVIntervalTree<String> initTree() {
        final SVIntervalTree<String> tree = new SVIntervalTree<>();
        final String[] genes = {"A", "B", "C", "D", "E"};
        Assert.assertEquals(transcriptionStartSites.length, genes.length);
        for ( int idx = 0; idx < genes.length; ++idx ) {
            tree.put(transcriptionStartSites[idx], genes[idx]);
        }
        return tree;
    }

    @DataProvider(name="nearestTSS")
    public Object[][] getNearestTSSData() {
        return new Object[][] {
                { new ClosedSVInterval(0, 1, 50), "A" },
                { new ClosedSVInterval(1, 105, 110), "B" },
                { new ClosedSVInterval(1, 155, 160), "B" },
                { new ClosedSVInterval(1, 160, 195), "C" },
                { new ClosedSVInterval(1, 3000, 4000), "D" },
                { new ClosedSVInterval(2, 33, 33), "E" },
                { new ClosedSVInterval(3, 900, 4000), null },
        };
    }

    @Test(dataProvider = "nearestTSS")
    public void testAnnotateNearestTranscriptionStartSite(
            final ClosedSVInterval variantInterval,
            final String expectedNearestTSSGene)
    {
        final String[] contigNames = {"chr1", "chr2", "chr3", "chr4"};
        final SVIntervalTree<String> transcriptionStartSiteTree = initTree();
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        SVAnnotate.annotateNearestTranscriptionStartSite(variantInterval.toSimpleInterval(contigNames),
                variantConsequenceDict, transcriptionStartSiteTree, 5000, variantInterval.getContig());
        if (isNull(expectedNearestTSSGene)) {
            Assert.assertNull(variantConsequenceDict.get(GATKSVVCFConstants.NEAREST_TSS));
        } else {
            Assert.assertEquals(variantConsequenceDict.get(GATKSVVCFConstants.NEAREST_TSS),
                    new HashSet<>(Arrays.asList(expectedNearestTSSGene)));
        }
    }

    private FeatureDataSource<GencodeGtfGeneFeature> loadToyGTFSource() {
        final File toyGTFFile = new File(getToolTestDataDir() + "unittest.gtf");
        return new FeatureDataSource<>(toyGTFFile);
    }

    @Test
    public void testGetTranscriptionStartSiteAndPromoterInterval() {
        Map<String, Integer> tssByGene = new HashMap<>();
        tssByGene.put("EMMA1", 100);
        tssByGene.put("EMMA2", 3000);

        Map<String, ClosedSVInterval> promoterByGene = new HashMap<>();
        promoterByGene.put("EMMA1", new ClosedSVInterval(0, 1,99));
        promoterByGene.put("EMMA2", new ClosedSVInterval(0, 3001, 4000));

        int promoterWindow = 1000;
        int contigID = 0;

        FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource();
        for (final GencodeGtfGeneFeature gene : toyGTFSource) {
            GencodeGtfTranscriptFeature transcript = gene.getTranscripts().get(0); // each gene only has one transcript
            String geneName = transcript.getGeneName();
            if (tssByGene.containsKey(geneName)) {
                int expectedTSS = tssByGene.get(geneName);
                Assert.assertEquals(expectedTSS, SVAnnotate.getTranscriptionStartSite(transcript));

                Assert.assertEquals(promoterByGene.get(geneName),
                        SVAnnotate.getPromoterInterval(transcript, promoterWindow, contigID));
            }
        }
    }


    private final Set<String> MSVExonOverlapClassifications = Sets.newHashSet(GATKSVVCFConstants.LOF,
            GATKSVVCFConstants.INT_EXON_DUP,
            GATKSVVCFConstants.DUP_PARTIAL,
            GATKSVVCFConstants.PARTIAL_EXON_DUP,
            GATKSVVCFConstants.COPY_GAIN,
            GATKSVVCFConstants.TSS_DUP);

    private GencodeGtfTranscriptFeature loadToyGtfTranscript() {
        FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource();
        // get only first gene, EMMA1, which has only one transcript
        final GencodeGtfGeneFeature toyGene = toyGTFSource.iterator().next();
        return toyGene.getTranscripts().get(0);
    }

    @DataProvider(name = "toyIntervalVariants")
    public Object[][] getToyIntervalVariantTestData() {
        return new Object[][] {
                { new SimpleInterval("chr1", 350, 550), GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 320, 360), GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 550, 750), GATKSVVCFConstants.LOF, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 450, 650), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 250, 850), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 150, 450), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 250, 350), GATKSVVCFConstants.PARTIAL_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 350, 650), GATKSVVCFConstants.PARTIAL_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 150, 250), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 650, 850), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 150, 750), GATKSVVCFConstants.INT_EXON_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 220, 280), GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC },
                { new SimpleInterval("chr1", 750, 950), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 50, 150), GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 50, 250), GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 50, 550), GATKSVVCFConstants.TSS_DUP, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 950, 1500), GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.UTR, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 650, 950), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 50, 2000), GATKSVVCFConstants.COPY_GAIN, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.INV_SPAN },
                { new SimpleInterval("chr1", 750, 790), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR },
                { new SimpleInterval("chr1", 350, 1200), GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF},
                { new SimpleInterval("chr1", 450, 1300), GATKSVVCFConstants.DUP_PARTIAL, GATKSVVCFConstants.LOF, GATKSVVCFConstants.MSV_EXON_OVERLAP, GATKSVVCFConstants.LOF}
        };
    }

    @Test(dataProvider = "toyIntervalVariants")
    public void testAnnotateIntervalSVTypes(
            final SimpleInterval toyVariant,
            final String expectedDuplicationConsequence,
            final String expectedDeletionConsequence,
            final String expectedCopyNumberVariantConsequence,
            final String expectedInversionConsequence
    ) {
        final GencodeGtfTranscriptFeature toyTranscript = loadToyGtfTranscript();

        final String actualDuplicationConsequence = SVAnnotate.annotateDuplication(toyVariant, toyTranscript);
        Assert.assertEquals(expectedDuplicationConsequence, actualDuplicationConsequence);

        final String actualDeletionConsequence = SVAnnotate.annotateDeletion(toyVariant, toyTranscript);
        Assert.assertEquals(expectedDeletionConsequence, actualDeletionConsequence);

        final String actualCopyNumberVariantConsequence = SVAnnotate.annotateCopyNumberVariant(toyVariant, toyTranscript, MSVExonOverlapClassifications);
        Assert.assertEquals(expectedCopyNumberVariantConsequence, actualCopyNumberVariantConsequence);

        final String actualInversionConsequence = SVAnnotate.annotateInversion(toyVariant, toyTranscript);
        Assert.assertEquals(expectedInversionConsequence, actualInversionConsequence);
    }

    @DataProvider(name = "toyPointVariants")
    public Object[][] getToyPointVariantTestData() {
        return new Object[][] {
                { new SimpleInterval("chr1", 150, 151), new SimpleInterval("chr1", 150, 150), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 250, 251), new SimpleInterval("chr1", 250, 250), GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.INTRONIC, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 350, 351), new SimpleInterval("chr1", 350, 350), GATKSVVCFConstants.LOF, GATKSVVCFConstants.BREAKEND_EXON, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 600, 601), new SimpleInterval("chr1", 600, 600), GATKSVVCFConstants.LOF, GATKSVVCFConstants.BREAKEND_EXON, GATKSVVCFConstants.LOF },
                { new SimpleInterval("chr1", 100, 101), new SimpleInterval("chr1", 100, 100), GATKSVVCFConstants.UTR, GATKSVVCFConstants.UTR, GATKSVVCFConstants.LOF }
        };
    }

    @Test(dataProvider = "toyPointVariants")
    public void testAnnotatePointSVTypes(
            final SimpleInterval toyTwoBaseVariant,
            final SimpleInterval toyPointVariant,
            final String expectedInsertionConsequence,
            final String expectedBreakendConsequence,
            final String expectedTranslocationVariantConsequence
    ) {
        final GencodeGtfTranscriptFeature toyTranscript = loadToyGtfTranscript();

        final String actualInsertionConsequence = SVAnnotate.annotateInsertion(toyTwoBaseVariant, toyTranscript);
        Assert.assertEquals(expectedInsertionConsequence, actualInsertionConsequence);

        // BND and CTX are annotated one breakpoint at a time
        final String actualBreakendConsequence = SVAnnotate.annotateBreakend(toyPointVariant, toyTranscript);
        Assert.assertEquals(expectedBreakendConsequence, actualBreakendConsequence);

        final String actualPointTranslocationConsequence = SVAnnotate.annotateTranslocation(toyPointVariant, toyTranscript);
        Assert.assertEquals(expectedTranslocationVariantConsequence, actualPointTranslocationConsequence);

        final String actualTwoBaseTranslocationConsequence = SVAnnotate.annotateTranslocation(toyTwoBaseVariant, toyTranscript);
        Assert.assertEquals(expectedTranslocationVariantConsequence, actualTwoBaseTranslocationConsequence);
    }

    @DataProvider(name = "toyComplexVariants")
    public Object[][] getToyComplexVariantTestData() {
        return new Object[][] {
                { "DUP_chr1:280-420", Sets.newHashSet(GATKSVVCFConstants.INT_EXON_DUP) },
                { "INV_chr1:90-1010,DUP_chr1:890-1010", Sets.newHashSet(GATKSVVCFConstants.INV_SPAN, GATKSVVCFConstants.DUP_PARTIAL) },
                { "DEL_chr1:250-450,INV_chr1:450-650,DUP_chr1:610-650", Sets.newHashSet(GATKSVVCFConstants.LOF, GATKSVVCFConstants.INTRONIC) }
        };
    }

    @Test(dataProvider = "toyComplexVariants")
    public void testAnnotateComplexEvents(
            final String cpxIntervalsString,
            final Set<String> expectedConsequences
    ) {
        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource();
        final Map<String, Set<String>> variantConsequenceDict = new HashMap<>();
        final Map<String,Integer> contigNameToID = createContigNameToIDMap(Arrays.asList("chr1"), Arrays.asList(0));
        final int promoterWindow = 1000;

        final SVAnnotate.GTFIntervalTreesContainer gtfTrees = SVAnnotate.buildIntervalTreesFromGTF(toyGTFSource, contigNameToID, promoterWindow);
        final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree = gtfTrees.gtfIntervalTree;

        final String[] cpxIntervalStrings = cpxIntervalsString.split(",");
        for (String cpxIntervalString : cpxIntervalStrings) {
            SVAnnotate.SVSegment cpxSegment = SVAnnotate.parseCPXIntervalString(cpxIntervalString);
            SVAnnotate.annotateGeneOverlaps(cpxSegment.interval, cpxSegment.intervalSVType,
                    variantConsequenceDict, MSVExonOverlapClassifications, contigNameToID, gtfIntervalTree);
        }

        Assert.assertEquals(expectedConsequences, variantConsequenceDict.keySet());
    }


    @Test
    public void testSortVariantConsequenceDict() {
        Map<String, Set<String>> before = new HashMap<>();
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "NOC2L");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "KLHL17");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "PLEKHN1");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "PERM1");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.DUP_PARTIAL, "SAMD11");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "HES4");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.TSS_DUP, "ISG15");

        Map<String, Object> expectedAfter = new HashMap<>();
        expectedAfter.put(GATKSVVCFConstants.DUP_PARTIAL, Arrays.asList("SAMD11"));
        expectedAfter.put(GATKSVVCFConstants.TSS_DUP, Arrays.asList("ISG15"));
        expectedAfter.put(GATKSVVCFConstants.LOF, Arrays.asList("HES4", "KLHL17", "NOC2L", "PERM1", "PLEKHN1"));

        Assert.assertEquals(expectedAfter, SVAnnotate.sortVariantConsequenceDict(before));
    }

    // create list of SV segments with SAME SVTYPE - convenience function for testing getSVSegments
    private List<SVAnnotate.SVSegment> createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType svType, SimpleInterval[] intervals) {
        List<SVAnnotate.SVSegment> segments = new ArrayList<>();
        for (SimpleInterval interval : intervals) {
            segments.add(new SVAnnotate.SVSegment(svType, interval));
        }
        return segments;
    }

    private void assertSegmentListEqual(List<SVAnnotate.SVSegment> segmentsA, List<SVAnnotate.SVSegment> segmentsB) {
        int lengthA = segmentsA.size();
        if (lengthA != segmentsB.size()) {
            Assert.fail("Segment lists differ in length");
        }
        for (int i = 0; i < lengthA; i++) {
            SVAnnotate.SVSegment segmentA = segmentsA.get(i);
            SVAnnotate.SVSegment segmentB = segmentsB.get(i);
            if (!segmentA.equals(segmentB)) {
                Assert.fail("Segment items differ");
            }
        }
    }

    private VariantContext createVariantContext(String chrom, int pos, int end, String chr2, Integer end2, String ref,
                                                String alt, Integer svLen, String strands, String cpxType,
                                                List<String> cpxIntervals) {
        final Map<String, Object> attributes = new HashMap<>();
        if (chr2 != null) {
            attributes.put(GATKSVVCFConstants.CONTIG2_ATTRIBUTE, chr2);
        }
        if (end2 != null) {
            attributes.put(GATKSVVCFConstants.END2_ATTRIBUTE, end2);
        }
        if (svLen != null) {
            attributes.put(GATKSVVCFConstants.SVLEN, svLen);
        }
        if (strands != null) {
            attributes.put(GATKSVVCFConstants.STRANDS_ATTRIBUTE, strands);
        }
        if (cpxType != null) {
            attributes.put(GATKSVVCFConstants.CPX_TYPE, cpxType);
        }
        if (cpxIntervals != null) {
            attributes.put(GATKSVVCFConstants.CPX_INTERVALS, cpxIntervals);
        }
        return new VariantContextBuilder()
                .source("source")
                .id("id")
                .chr(chrom)
                .start(pos)
                .stop(end)
                .alleles(Arrays.asList(ref != null ? Allele.create(ref, true) : Allele.REF_N,
                        alt != null ? Allele.create(alt, false) : Allele.ALT_N))
                .attributes(attributes)
                .make();
    }

    @DataProvider(name = "typesAndSegments")
    public Object[][] getSVTypesAndSegmentsTestData() {
        return new Object[][] {
                { createVariantContext("chr2", 86263976, 86263977, "chr19", 424309, "N",
                        "<CTX>", null, null, "CTX_PP/QQ", null),
                        SVAnnotate.StructuralVariantAnnotationType.CTX,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval[]{ new SimpleInterval("chr2", 86263976, 86263977),
                                        new SimpleInterval("chr19", 424309, 424310) }),
                        null },
                { createVariantContext("chr2", 86263976, 86263976, null, 424309, "G",
                        "G]chr19:424309]", null, null,"CTX_PP/QQ", null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr2", 86263976, 86263976))),
                        null},
                { createVariantContext("chr2", 86263977, 86263977, null, 424309, "A",
                        "[chr19:424310[A", null, null, "CTX_PP/QQ", null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr2", 86263977, 86263977))),
                        null },
                { createVariantContext("chr2", 205522308, 205522384, "chr2", null, "N",
                        "<INV>", 76, null, null, null),
                        SVAnnotate.StructuralVariantAnnotationType.INV,
                        Arrays.asList(new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.INV,
                                new SimpleInterval("chr2", 205522308, 205522384))),
                        null },
                { createVariantContext("chr19", 424309, 424309, null, 424309, "T",
                        "T]chr2:86263976]", null, null, "CTX_PP/QQ", null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr19", 424309, 424309))),
                        null },
                { createVariantContext("chr19", 424310, 424310, null, 424309, "C",
                        "[chr2:86263977[C", null, null, "CTX_PP/QQ", null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        Arrays.asList(new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.CTX,
                                new SimpleInterval("chr19", 424310, 424310))),
                        null },
                { createVariantContext("chr22", 10510000, 10694100, "chr22", null, "N",
                        "<DEL>", 184100, null, null, null),
                        SVAnnotate.StructuralVariantAnnotationType.DEL,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DEL,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10510000, 10694100)}),
                        null},
                { createVariantContext("chr22", 10524000, 10710000, "chr22", null, "N",
                        "<DUP>", 186000, null, null, null),
                        SVAnnotate.StructuralVariantAnnotationType.DUP,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DUP,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10524000, 10710000)}),
                        null },
                { createVariantContext("chr22", 10532563, 10532611, "chr22", null, "N",
                        "<INS:ME:ALU>", 245, null, null, null),
                        SVAnnotate.StructuralVariantAnnotationType.INS,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.INS,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10532563, 10532564)}),
                        null },
                { createVariantContext("chr22", 10572758, 10572788, "chr22", null, "N",
                        "<INS>", 57, null, null, null),
                        SVAnnotate.StructuralVariantAnnotationType.INS,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.INS,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10572758, 10572759)}),
                        null },
                { createVariantContext("chr22", 10717890, 10717890, "chr22", null, "N",
                        "<BND>", 5170, "-+", null, null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10717890, 10717890),
                                        new SimpleInterval("chr22", 10723060, 10723060) }),
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DUP,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10717890, 10723060)}) },
                { createVariantContext("chr22", 10774600, 10784500, "chr22", null, "N",
                        "<CNV>", 9900, null, null, null),
                        SVAnnotate.StructuralVariantAnnotationType.CNV,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CNV,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10774600, 10784500)}),
                        null },
                { createVariantContext("chr22", 10930458, 10930458, "chr22", 11564561, "N",
                        "<BND>", 634103, "--", null, null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 10930458, 10930458),
                                        new SimpleInterval("chr22", 11564561, 11564561) }),
                        null },
                { createVariantContext("chr22", 17636024, 17636024, "chr22", null, "N",
                        "<BND>", 10709, "+-", null, null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 17636024, 17636024),
                                        new SimpleInterval("chr22", 17646733, 17646733) }),
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DEL,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 17636024, 17646733)}) },
                { createVariantContext("chr22", 18971159, 18971435, "chr22", null, "N",
                        "<CPX>", 386, null, "dDUP", Arrays.asList("INV_chr22:20267228-20267614","DUP_chr22:20267228-20267614")),
                        SVAnnotate.StructuralVariantAnnotationType.CPX,
                        Arrays.asList(
                                new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.INV,
                                        new SimpleInterval("chr22", 20267228, 20267614)),
                                new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.DUP,
                                        new SimpleInterval("chr22", 20267228, 20267614)),
                                new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.INS,
                                        new SimpleInterval("chr22", 18971159, 18971160))),
                        null },
                { createVariantContext("chr22", 22120897, 22120897, "chrX", 126356858, "N",
                        "<BND>", -1, "++", null, null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22120897, 22120897),
                                        new SimpleInterval("chrX", 126356858, 126356858) }),
                        null },
                { createVariantContext("chr22", 22196261, 22196261, "chr22", null, "N",
                        "<BND>", 708725, "+-", null, null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22196261, 22196261),
                                        new SimpleInterval("chr22", 22904986, 22904986) }),
                        null },
                { createVariantContext("chr22", 22196261, 22196261, null, null, "A",
                        "A[chr22:22904986[", null, "+-", null, null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22196261, 22196261) }),
                        null },
                { createVariantContext("chr22", 22904986, 22904986, null, null, "T",
                        "]chr22:22196261]T", null, "+-", null, null),
                        SVAnnotate.StructuralVariantAnnotationType.BND,
                        createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                                new SimpleInterval[]{ new SimpleInterval("chr22", 22904986, 22904986) }),
                        null },
                { createVariantContext("chr22", 36533058, 36538234, "chr22", null, "N",
                        "<CPX>", 5176, null, "dupINV", Arrays.asList("DUP_chr22:36533058-36533299","INV_chr22:36533058-36538234")),
                        SVAnnotate.StructuralVariantAnnotationType.CPX,
                        Arrays.asList(
                                new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.DUP,
                                        new SimpleInterval("chr22", 36533058, 36533299)),
                                new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.INV,
                                        new SimpleInterval("chr22", 36533058, 36538234))),
                        null }
        };
    }

    @Test(dataProvider = "typesAndSegments")
    public void testGetSVTypeAndSegments(
            final VariantContext variant,
            final SVAnnotate.StructuralVariantAnnotationType expectedSVType,
            final List<SVAnnotate.SVSegment> expectedSVSegments,
            final List<SVAnnotate.SVSegment> expectedSVSegmentsWithBNDOverlap
    ) {
        SVAnnotate.StructuralVariantAnnotationType actualSVType = SVAnnotate.getSVType(variant);
        Assert.assertEquals(actualSVType, expectedSVType);

        List<SVAnnotate.SVSegment> actualSegments = SVAnnotate.getSVSegments(variant, actualSVType, -1);
        assertSegmentListEqual(actualSegments, expectedSVSegments);

        List<SVAnnotate.SVSegment> actualSegmentsWithBNDOverlap = SVAnnotate.getSVSegments(variant, actualSVType,
                15000);
        assertSegmentListEqual(actualSegmentsWithBNDOverlap,
                expectedSVSegmentsWithBNDOverlap != null ? expectedSVSegmentsWithBNDOverlap : expectedSVSegments);
    }

    private Map<String, Object> createAttributesMap(List<String> keys, List<Object> values) {
        Map<String, Object> attributes = new HashMap<>();
        int len = keys.size();
        if (len != values.size()) {
            throw new IllegalArgumentException("Length of keys list != length of values list");
        }
        for (int i = 0; i < len; i++) {
            Object value = values.get(i);
            attributes.put(keys.get(i), value.getClass() == String.class ? Arrays.asList(value) : value);
        }
        return attributes;
    }

    @DataProvider(name = "toyIntegratedVariants")
    public Object[][] getAnnotateStructuralVariantTestData() {
        return new Object[][]{
                { createVariantContext("chr1", 10, 50, null, null, null,
                        "<CNV>", 40, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PROMOTER, GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", true)) },
                { createVariantContext("chr1", 50, 450, null, null, null,
                        "<DEL>", 400, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.LOF, GATKSVVCFConstants.NONCODING_BREAKPOINT,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA1", "DNase", false)) },
                { createVariantContext("chr1", 1100, 1700, null, null, null,
                        "<INV>", 600, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_SPAN, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("Enhancer", "EMMA1", true)) },
                { createVariantContext("chr1", 1300, 1800, null, null, null,
                        "<INS:LINE1>", 42, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.NONCODING_BREAKPOINT, GATKSVVCFConstants.NEAREST_TSS,
                                        GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("Enhancer", "EMMA1", true)) },
                { createVariantContext("chr1", 2000, 2200, null, null, null,
                        "<DUP>", 42, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.DUP_PARTIAL,GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA2", false)) },
                { createVariantContext("chr1", 3100, 3300, null, null, null,
                        "<DUP>", 200, null, null, null),
                        createAttributesMap(
                                Arrays.asList(GATKSVVCFConstants.PROMOTER,GATKSVVCFConstants.INTERGENIC),
                                Arrays.asList("EMMA2", true)) }
        };
    }

    private FeatureDataSource<FullBEDFeature> loadTinyNoncodingBEDSource() {
        return new FeatureDataSource<>(getToolTestDataDir() + "noncoding.unittest.bed");
    }

    /**
     * Tests full integrated annotateStructuralVariant process, with a focus on noncoding & promoter annotation
     * and logic controlling when to annotate TSS and intergenic, since protein-coding and TSS annotation
     * is tested elsewhere
     * @param variant
     * @param expectedAttributes
     */
    @Test(dataProvider = "toyIntegratedVariants")
    public void testAnnotateStructuralVariant(
            final VariantContext variant,
            final Map<String, Object> expectedAttributes
    ) {
        final Map<String,Integer> contigNameToID = createContigNameToIDMap(Arrays.asList("chr1"), Arrays.asList(0));
        final String[] contigIDToName = new String[]{"chr1"};
        final int promoterWindow = 200;
        final int maxBreakendLen = -1;
        final int maxContigLength = 5000;

        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource();
        final SVAnnotate.GTFIntervalTreesContainer gtfTrees = SVAnnotate.buildIntervalTreesFromGTF(toyGTFSource, contigNameToID, promoterWindow);
        final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree = gtfTrees.gtfIntervalTree;
        final SVIntervalTree<String> promoterIntervalTree = gtfTrees.promoterIntervalTree;
        final SVIntervalTree<String> transcriptionStartSiteTree = gtfTrees.transcriptionStartSiteTree;
        final FeatureDataSource<FullBEDFeature> tinyNoncodingBedSource = loadTinyNoncodingBEDSource();
        final SVIntervalTree<String> nonCodingIntervalTree =
                SVAnnotate.buildIntervalTreeFromBED(tinyNoncodingBedSource, contigNameToID);

        final Map<String, Object> actualAttributes = SVAnnotate.annotateStructuralVariant(variant, gtfIntervalTree, promoterIntervalTree,
                transcriptionStartSiteTree, nonCodingIntervalTree, MSVExonOverlapClassifications, contigNameToID,
                contigIDToName, maxBreakendLen, maxContigLength);

        Assert.assertEquals(expectedAttributes, actualAttributes);
    }

    private Map<String, Integer> createContigNameToIDMap(List<String> names, List<Integer> IDs) {
        final Map<String,Integer> contigNameToID = new HashMap<>();
        for (int i = 0 ; i < names.size() ; i++) {
            contigNameToID.put(names.get(i), IDs.get(i));
        }
        return contigNameToID;
    }

    @DataProvider(name = "buildIntervalTrees")
    public Object[][] getBuildIntervalTreesTestData() {
        return new Object[][] {
                { createContigNameToIDMap(Arrays.asList("chr1"), Arrays.asList(0)), 2, 2 },
                { createContigNameToIDMap(Arrays.asList("chr1", "chr2"), Arrays.asList(0, 1)), 3, 3 }
        };
    }

    @Test(dataProvider = "buildIntervalTrees")
    public void testIgnoreUnknownContigsWhenBuildingIntervalTrees(
            final Map<String, Integer> contigNameToID,
            final int expectedBEDTreeSize,
            final int expectedTranscriptTreeSize
    ) {
        final FeatureDataSource<FullBEDFeature> tinyNoncodingBedSource = loadTinyNoncodingBEDSource();
        final SVIntervalTree<String> nonCodingIntervalTree =
                SVAnnotate.buildIntervalTreeFromBED(tinyNoncodingBedSource, contigNameToID);
        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource();
        final SVAnnotate.GTFIntervalTreesContainer gtfTrees = SVAnnotate.buildIntervalTreesFromGTF(toyGTFSource, contigNameToID, 100);
        final SVIntervalTree<GencodeGtfTranscriptFeature> gtfIntervalTree = gtfTrees.gtfIntervalTree;
        // check size to ensure contigs not included in the map are excluded from the interval tree successfully
        Assert.assertEquals(nonCodingIntervalTree.size(), expectedBEDTreeSize);
        Assert.assertEquals(gtfIntervalTree.size(), expectedTranscriptTreeSize);

    }

    // interval trees? getContigIDFromName?
    // variantOverlapsTranscriptionStartSite
}