package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.variant.variantcontext.VariantContext;
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
import java.sql.Array;
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
    private static ClosedSVInterval[] transcriptionStartSites = {
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

    private final FeatureDataSource<GencodeGtfGeneFeature> loadToyGTFSource() {
        final File toyGTFFile = new File(getToolTestDataDir() + "EMMA1.gtf");
        final FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = new FeatureDataSource<>(toyGTFFile);
        return toyGTFSource;
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

    private final GencodeGtfTranscriptFeature loadToyGtfTranscript() {
        FeatureDataSource<GencodeGtfGeneFeature> toyGTFSource = loadToyGTFSource();
        // get only first gene, EMMA1, which has only one transcript
        final GencodeGtfGeneFeature toyGene = toyGTFSource.iterator().next();
        final GencodeGtfTranscriptFeature toyTranscript = toyGene.getTranscripts().get(0);
        return toyTranscript;
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
        final Map<String,Integer> contigNameToID = new HashMap<>();
        contigNameToID.put("chr1", 0);
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
    public void testFormatVariantConsequenceDict() {
        Map<String, Set<String>> before = new HashMap<>();
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "NOC2L");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "KLHL17");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "PLEKHN1");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "PERM1");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.DUP_PARTIAL, "SAMD11");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.LOF, "HES4");
        SVAnnotate.updateVariantConsequenceDict(before, GATKSVVCFConstants.TSS_DUP, "ISG15");

        Map<String, String> expectedAfter = new HashMap<>();
        expectedAfter.put(GATKSVVCFConstants.DUP_PARTIAL, "SAMD11");
        expectedAfter.put(GATKSVVCFConstants.TSS_DUP, "ISG15");
        expectedAfter.put(GATKSVVCFConstants.LOF, "HES4,KLHL17,NOC2L,PERM1,PLEKHN1");

        Assert.assertEquals(expectedAfter, SVAnnotate.formatVariantConsequenceDict(before));
    }

    // create list of SV segments with SAME SVTYPE - convenience function for testing
    private List<SVAnnotate.SVSegment> createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType svType, SimpleInterval[] intervals) {
        List<SVAnnotate.SVSegment> segments = new ArrayList<>();
        for (SimpleInterval interval : intervals) {
            segments.add(new SVAnnotate.SVSegment(svType, interval));
        }
        return segments;
    }

    private Map<String, SVAnnotate.StructuralVariantAnnotationType> createSVTypeByVariantIDMap() {
        Map<String, SVAnnotate.StructuralVariantAnnotationType> SVTypeByVariantID = new HashMap<>();
        SVTypeByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1", SVAnnotate.StructuralVariantAnnotationType.CTX);
        SVTypeByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M1", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M3", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_INV_chr2_5", SVAnnotate.StructuralVariantAnnotationType.INV);
        SVTypeByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M2", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M4", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_DEL_chr22_1", SVAnnotate.StructuralVariantAnnotationType.DEL);
        SVTypeByVariantID.put("ref_panel_1kg_v1_DUP_chr22_1", SVAnnotate.StructuralVariantAnnotationType.DUP);
        SVTypeByVariantID.put("ref_panel_1kg_v1_INS_chr22_1", SVAnnotate.StructuralVariantAnnotationType.INS);
        SVTypeByVariantID.put("ref_panel_1kg_v1_INS_chr22_2", SVAnnotate.StructuralVariantAnnotationType.INS);
        SVTypeByVariantID.put("ref_panel_1kg_v1_BND_chr22_1", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_CNV_chr22_1", SVAnnotate.StructuralVariantAnnotationType.CNV);
        SVTypeByVariantID.put("ref_panel_1kg_v1_BND_chr22_3", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_BND_chr22_7", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_CPX_chr22_1", SVAnnotate.StructuralVariantAnnotationType.CPX);
        SVTypeByVariantID.put("ref_panel_1kg_v1_BND_chr22_14", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_BND_chr22_16", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_BND_chr22_16_M1", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_BND_chr22_16_M2", SVAnnotate.StructuralVariantAnnotationType.BND);
        SVTypeByVariantID.put("ref_panel_1kg_v1_CPX_chr22_3", SVAnnotate.StructuralVariantAnnotationType.CPX);
        return SVTypeByVariantID;
    }

    private Map<String, List<SVAnnotate.SVSegment>> createSegmentsByVariantIDMap() {
        Map<String, List<SVAnnotate.SVSegment>> segmentsByVariantID = new HashMap<>();
        segmentsByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CTX,
                        new SimpleInterval[]{ new SimpleInterval("chr2", 86263976, 86263977),
                                new SimpleInterval("chr19", 424309, 424310) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CTX,
                        new SimpleInterval[]{ new SimpleInterval("chr2", 86263976, 86263976) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M3",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CTX,
                        new SimpleInterval[]{ new SimpleInterval("chr2", 86263977, 86263977) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_INV_chr2_5",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.INV,
                        new SimpleInterval[]{ new SimpleInterval("chr2", 205522308, 205522384) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M2",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CTX,
                        new SimpleInterval[]{ new SimpleInterval("chr19", 424309, 424309) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_CTX_chr2_1_M4",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CTX,
                        new SimpleInterval[]{ new SimpleInterval("chr19", 424310, 424310) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_DEL_chr22_1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DEL,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10510000, 10694100) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_DUP_chr22_1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DUP,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10524000, 10710000) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_INS_chr22_1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.INS,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10532563, 10532564) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_INS_chr22_2",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.INS,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10572758, 10572759) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_BND_chr22_1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10717890, 10717890),
                                new SimpleInterval("chr22", 10723060, 10723060) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_CNV_chr22_1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.CNV,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10774600, 10784500) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_BND_chr22_3",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10930458, 10930458),
                                new SimpleInterval("chr22", 11564561, 11564561) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_BND_chr22_7",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 17636024, 17636024),
                                new SimpleInterval("chr22", 17646733, 17646733) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_CPX_chr22_1",
                Arrays.asList(new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.INV,
                                new SimpleInterval("chr22", 20267228, 20267614)),
                        new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.DUP,
                                new SimpleInterval("chr22", 20267228, 20267614)),
                        new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.INS,
                                new SimpleInterval("chr22", 18971159, 18971160))));
        segmentsByVariantID.put("ref_panel_1kg_v1_BND_chr22_14",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 22120897, 22120897),
                                new SimpleInterval("chrX", 126356858, 126356858) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_BND_chr22_16",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 22196261, 22196261),
                                new SimpleInterval("chr22", 22904986, 22904986) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_BND_chr22_16_M1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 22196261, 22196261) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_BND_chr22_16_M2",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.BND,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 22904986, 22904986) }));
        segmentsByVariantID.put("ref_panel_1kg_v1_CPX_chr22_3",
                Arrays.asList(new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.DUP,
                                new SimpleInterval("chr22", 36533058, 36533299)),
                        new SVAnnotate.SVSegment(SVAnnotate.StructuralVariantAnnotationType.INV,
                                new SimpleInterval("chr22", 36533058, 36538234))));
        return segmentsByVariantID;
    }

    // create copy of segments answers and update as needed if maxBreakendLenForOverlapAnnotation = 15000
    private Map<String, List<SVAnnotate.SVSegment>> createSegmentsByVariantIDMapForBNDOverlap(Map<String, List<SVAnnotate.SVSegment>> segmentsByVariantID) {
        Map<String, List<SVAnnotate.SVSegment>> updated = new HashMap<>(segmentsByVariantID);
        updated.put("ref_panel_1kg_v1_BND_chr22_1",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DUP,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 10717890, 10723060)}));
        updated.put("ref_panel_1kg_v1_BND_chr22_7",
                createListOfSVSegments(SVAnnotate.StructuralVariantAnnotationType.DEL,
                        new SimpleInterval[]{ new SimpleInterval("chr22", 17636024, 17646733)}));
        return updated;
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

    @Test
    public void testGetSVTypeAndSegments() {
        Map<String, SVAnnotate.StructuralVariantAnnotationType> SVTypeByVariantID = createSVTypeByVariantIDMap();
        Map<String, List<SVAnnotate.SVSegment>> segmentsByVariantID = createSegmentsByVariantIDMap();
        Map<String, List<SVAnnotate.SVSegment>> segmentsWithBNDOverlapsByVariantID = createSegmentsByVariantIDMapForBNDOverlap(segmentsByVariantID);
        final File typesVCFFile = new File(getToolTestDataDir() + "types.vcf");
        final FeatureDataSource<VariantContext> typesVariants = new FeatureDataSource<>(typesVCFFile);
        for (VariantContext variant : typesVariants) {
            final String variantID = variant.getID();
            SVAnnotate.StructuralVariantAnnotationType actualSVType = SVAnnotate.getSVType(variant);
            Assert.assertEquals(actualSVType, SVTypeByVariantID.get(variantID));

            List<SVAnnotate.SVSegment> actualSegments = SVAnnotate.getSVSegments(variant, actualSVType, -1);
            // Assert.assertEquals(actualSegments, segmentsByVariantID.get(variantID));
            List<SVAnnotate.SVSegment> expectedSegments = segmentsByVariantID.get(variantID);
            assertSegmentListEqual(actualSegments, expectedSegments);

            List<SVAnnotate.SVSegment> actualSegmentsWithBNDOverlap = SVAnnotate.getSVSegments(variant, actualSVType, 15000);
            // Assert.assertEquals(actualSegmentsWithBNDOverlap, segmentsByVariantID.get(variantID));
            assertSegmentListEqual(actualSegmentsWithBNDOverlap, segmentsWithBNDOverlapsByVariantID.get(variantID));
        }
    }

    // noncoding and promoter annotation unit tests
    // getSVSegments and getSVType
    // interval trees? getContigIDFromName?
    // variantOverlapsTranscriptionStartSite
}