package org.broadinstitute.hellbender.tools.spark.sv.discovery.readdepth;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class LargeSimpleSVFactoryTest extends BaseTest {

    private static LargeSimpleSVFactory getEmptyTandemDuplicationFactory(final List<EvidenceTargetLink> intrachromosomalEvidenceTargetLinks,
                                                                         final List<CopyRatio> copyRatioList,
                                                                         final List<CalledCopyRatioSegment> copyRatioSegmentList,
                                                                         final SAMSequenceDictionary dictionary) {
        final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree = new SVIntervalTree<>();
        for (final EvidenceTargetLink link : intrachromosomalEvidenceTargetLinks) {
            final SVInterval leftInterval = link.getPairedStrandedIntervals().getLeft().getInterval();
            final SVInterval rightInterval = link.getPairedStrandedIntervals().getRight().getInterval();
            intrachromosomalLinkTree.put(leftInterval, link);
            intrachromosomalLinkTree.put(rightInterval, link);
        }
        final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree = new SVIntervalTree<>();
        final SVIntervalTree<GATKRead> contigTree = new SVIntervalTree<>();
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection();

        final SampleLocatableMetadata sampleLocatableMetadata = new SimpleSampleLocatableMetadata("test", dictionary);
        final CalledCopyRatioSegmentCollection calledCopyRatioSegmentCollection = new CalledCopyRatioSegmentCollection(sampleLocatableMetadata, copyRatioSegmentList);
        final CopyRatioCollection copyRatioSegmentCollection = new CopyRatioCollection(sampleLocatableMetadata, copyRatioList);

        final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector = calledCopyRatioSegmentCollection.getOverlapDetector();
        final OverlapDetector<CopyRatio> copyRatioOverlapDetector = copyRatioSegmentCollection.getOverlapDetector();

        return new LargeTandemDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree, contigTree, arguments, copyRatioSegmentOverlapDetector, copyRatioOverlapDetector, dictionary);
    }

    @Test(groups = "sv")
    public void testCall() {

        //Test case with amplification at contig0:300-350
        final List<Double> rawCopyRatios = new ArrayList<>(100);
        for (int i = 0; i < 100; i++) {
            if (i >= 30 && i < 35) {
                rawCopyRatios.add(1.5);
            } else {
                rawCopyRatios.add(1.0);
            }
        }
        final int binSize = 10;
        final List<CopyRatio> copyRatiosList = getCopyRatioBins(rawCopyRatios, "contig0", 0, binSize);
        final List<CalledCopyRatioSegment> calledCopyRatioSegmentList = Arrays.asList(
                new CalledCopyRatioSegment(new CopyRatioSegment(new SimpleInterval("contig0", 1, 300), 30, 0.0), CalledCopyRatioSegment.Call.NEUTRAL),
                new CalledCopyRatioSegment(new CopyRatioSegment(new SimpleInterval("contig0", 301, 350), 5, 0.585), CalledCopyRatioSegment.Call.AMPLIFICATION),
                new CalledCopyRatioSegment(new CopyRatioSegment(new SimpleInterval("contig0", 351, binSize * rawCopyRatios.size()), 65, 0.0), CalledCopyRatioSegment.Call.NEUTRAL));
        final List<EvidenceTargetLink> intrachromosomalEvidenceTargetLinks = Arrays.asList(
                new EvidenceTargetLink(new StrandedInterval(new SVInterval(0, 290, 310), false),
                        new StrandedInterval(new SVInterval(0, 340, 360), true),
                        0, 5, new HashSet<>(Arrays.asList("read1", "read2", "read3", "read4", "read5")), Collections.emptySet()));
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("contig0", 2000));
        final LargeSimpleSVFactory factory = getEmptyTandemDuplicationFactory(intrachromosomalEvidenceTargetLinks, copyRatiosList, calledCopyRatioSegmentList, dictionary);

        final LargeSimpleSV emptyResult = factory.call(new SVInterval(0, 600, 700), new SVInterval(0, 1100, 1200), new SVInterval(0, 650, 1150), null, 0);
        Assert.assertNull(emptyResult);

        final LargeSimpleSV result = factory.call(new SVInterval(0, 285, 295), new SVInterval(0, 345, 355), new SVInterval(0, 290, 350), null, 10);
        Assert.assertEquals(result, new LargeSimpleSV(SimpleSVType.TYPES.DUP, 290, 350, 0, 5, 0, 0, 0, null));
    }

    @Test(groups = "sv")
    public void testCountSplitReads() {

        final List<EvidenceTargetLink> list = new ArrayList<>();
        Assert.assertEquals(LargeSimpleSVFactory.countSplitReads(list), 0);

        final StrandedInterval source = new StrandedInterval(new SVInterval(0, 100, 200), true);
        final StrandedInterval target = new StrandedInterval(new SVInterval(0, 300, 400), false);
        list.add(new EvidenceTargetLink(source, target, 0, 1, Collections.singleton("read1"), Collections.emptySet()));
        Assert.assertEquals(LargeSimpleSVFactory.countSplitReads(list), 0);

        list.add(new EvidenceTargetLink(source, target, 1, 0, Collections.emptySet(), Collections.singleton("read2")));
        Assert.assertEquals(LargeSimpleSVFactory.countSplitReads(list), 1);

        list.add(new EvidenceTargetLink(source, target, 1, 0, Collections.emptySet(), Collections.singleton("read3")));
        Assert.assertEquals(LargeSimpleSVFactory.countSplitReads(list), 2);

        list.add(new EvidenceTargetLink(source, target, 1, 0, Collections.emptySet(), Collections.singleton("read3")));
        Assert.assertEquals(LargeSimpleSVFactory.countSplitReads(list), 2);
    }

    @Test(groups = "sv")
    public void testReadPairEvidence() {

        final List<EvidenceTargetLink> list = new ArrayList<>();
        Assert.assertEquals(LargeSimpleSVFactory.countReadPairs(list), 0);

        final StrandedInterval source = new StrandedInterval(new SVInterval(0, 100, 200), true);
        final StrandedInterval target = new StrandedInterval(new SVInterval(0, 300, 400), false);
        list.add(new EvidenceTargetLink(source, target, 1, 0, Collections.emptySet(), Collections.singleton("read1")));
        Assert.assertEquals(LargeSimpleSVFactory.countReadPairs(list), 0);

        list.add(new EvidenceTargetLink(source, target, 0, 1, Collections.singleton("read2"), Collections.emptySet()));
        Assert.assertEquals(LargeSimpleSVFactory.countReadPairs(list), 1);

        list.add(new EvidenceTargetLink(source, target, 0, 1, Collections.singleton("read3"), Collections.emptySet()));
        Assert.assertEquals(LargeSimpleSVFactory.countReadPairs(list), 2);

        list.add(new EvidenceTargetLink(source, target, 0, 1, Collections.singleton("read3"), Collections.emptySet()));
        Assert.assertEquals(LargeSimpleSVFactory.countReadPairs(list), 2);
    }

    @Test(groups = "sv")
    public void testGetCopyRatiosOnInterval() {
        final List<CopyRatio> copyRatioList = new ArrayList<>();
        for (int i = 1; i < 1000; i += 10) {
            copyRatioList.add(new CopyRatio(new SimpleInterval("contig0", i, i + 9), 1.0));
        }

        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("contig0", 2000));
        final SampleLocatableMetadata sampleLocatableMetadata = new SimpleSampleLocatableMetadata("test", dictionary);
        final CopyRatioCollection copyRatioCollection = new CopyRatioCollection(sampleLocatableMetadata, copyRatioList);
        final OverlapDetector<CopyRatio> overlapDetector = copyRatioCollection.getOverlapDetector();

        final SVInterval interval1 = new SVInterval(0, 101, 200);
        final List<CopyRatio> result1 = LargeSimpleSVFactory.getCopyRatiosOnInterval(interval1, overlapDetector, 0, dictionary);
        final List<CopyRatio> expectedResult1 = copyRatioList.subList(10, 20);
        Assert.assertEquals(result1, expectedResult1);

        final SVInterval interval2 = new SVInterval(0, 1100, 1200);
        final List<CopyRatio> result2 = LargeSimpleSVFactory.getCopyRatiosOnInterval(interval2, overlapDetector, 0, dictionary);
        Assert.assertTrue(result2.isEmpty());

        final List<CopyRatio> result3 = LargeSimpleSVFactory.getCopyRatiosOnInterval(interval1, overlapDetector, 1, dictionary);
        final List<CopyRatio> expectedResult3 = copyRatioList.subList(11, 19);
        Assert.assertEquals(result3, expectedResult3);

        final List<CopyRatio> result4 = LargeSimpleSVFactory.getCopyRatiosOnInterval(interval1, overlapDetector, 5, dictionary);
        Assert.assertTrue(result4.isEmpty());

        final List<CopyRatio> result5 = LargeSimpleSVFactory.getCopyRatiosOnInterval(interval1, overlapDetector, 6, dictionary);
        Assert.assertTrue(result5.isEmpty());
    }

    @Test(groups = "sv")
    public void testCountUniqueCounterEvidence() {

        final List<EvidenceTargetLink> evidenceTargetLinks = new ArrayList<>();
        final List<EvidenceTargetLink> counterEvidenceTargetLinks = new ArrayList<>();
        Assert.assertEquals(LargeSimpleSVFactory.countUniqueCounterEvidence(counterEvidenceTargetLinks, evidenceTargetLinks, 1, EvidenceTargetLink::getReadPairTemplateNames), 0);

        final StrandedInterval source = new StrandedInterval(new SVInterval(0, 100, 200), true);
        final StrandedInterval target = new StrandedInterval(new SVInterval(0, 300, 400), false);
        counterEvidenceTargetLinks.add(new EvidenceTargetLink(source, target, 0, 1, Collections.singleton("read1"), Collections.emptySet()));
        Assert.assertEquals(LargeSimpleSVFactory.countUniqueCounterEvidence(counterEvidenceTargetLinks, evidenceTargetLinks, 1, EvidenceTargetLink::getReadPairTemplateNames), 1);
        Assert.assertEquals(LargeSimpleSVFactory.countUniqueCounterEvidence(counterEvidenceTargetLinks, evidenceTargetLinks, 1, EvidenceTargetLink::getSplitReadTemplateNames), 0);
        Assert.assertEquals(LargeSimpleSVFactory.countUniqueCounterEvidence(counterEvidenceTargetLinks, evidenceTargetLinks, 2, EvidenceTargetLink::getReadPairTemplateNames), 0);

        evidenceTargetLinks.add(new EvidenceTargetLink(source, target, 1, 0, Collections.emptySet(), Collections.singleton("read1")));
        Assert.assertEquals(LargeSimpleSVFactory.countUniqueCounterEvidence(counterEvidenceTargetLinks, evidenceTargetLinks, 1, EvidenceTargetLink::getReadPairTemplateNames), 0);

        evidenceTargetLinks.remove(evidenceTargetLinks.size() - 1);
        evidenceTargetLinks.add(new EvidenceTargetLink(source, target, 0, 1, Collections.singleton("read1"), Collections.emptySet()));
        Assert.assertEquals(LargeSimpleSVFactory.countUniqueCounterEvidence(counterEvidenceTargetLinks, evidenceTargetLinks, 1, EvidenceTargetLink::getReadPairTemplateNames), 0);

        counterEvidenceTargetLinks.add(new EvidenceTargetLink(source, target, 2, 0, Collections.emptySet(), new HashSet<>(Arrays.asList("read2", "read3"))));
        Assert.assertEquals(LargeSimpleSVFactory.countUniqueCounterEvidence(counterEvidenceTargetLinks, evidenceTargetLinks, 1, EvidenceTargetLink::getSplitReadTemplateNames), 2);
    }

    @Test(groups = "sv")
    public void testGetMatchingLinks() {
        final SVIntervalTree<EvidenceTargetLink> tree = new SVIntervalTree<>();

        final EvidenceTargetLink link = new EvidenceTargetLink(
                new StrandedInterval(new SVInterval(0, 100, 200), true),
                new StrandedInterval(new SVInterval(0, 600, 700), false),
                0, 0, Collections.emptySet(), Collections.emptySet());
        tree.put(new SVInterval(0, 100, 200), link);
        tree.put(new SVInterval(0, 600, 700), link);

        Assert.assertTrue(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(0, 10, 50),
                new SVInterval(0, 250, 300), tree).isEmpty());
        Assert.assertTrue(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(1, 100, 200),
                new SVInterval(0, 600, 700), tree).isEmpty());
        Assert.assertTrue(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(0, 90, 190),
                new SVInterval(0, 800, 900), tree).isEmpty());
        Assert.assertTrue(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(0, 300, 400),
                new SVInterval(0, 600, 700), tree).isEmpty());

        Assert.assertEquals(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(0, 150, 250),
                new SVInterval(0, 550, 650), tree), Collections.singletonList(link));
        Assert.assertEquals(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(0, 150, 200),
                new SVInterval(0, 550, 650), tree), Collections.singletonList(link));
        Assert.assertEquals(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(0, 199, 250),
                new SVInterval(0, 550, 650), tree), Collections.singletonList(link));

        final EvidenceTargetLink link2 = new EvidenceTargetLink(
                new StrandedInterval(new SVInterval(0, 150, 250), true),
                new StrandedInterval(new SVInterval(0, 650, 750), false),
                0, 0, Collections.emptySet(), Collections.emptySet());
        tree.put(new SVInterval(0, 150, 250), link2);
        tree.put(new SVInterval(0, 650, 750), link2);
        Assert.assertEquals(LargeSimpleSVFactory.getMatchingLinks(new SVInterval(0, 150, 250),
                new SVInterval(0, 600, 700), tree), Arrays.asList(link, link2));
    }

    @Test(groups = "sv")
    public void testLocalOverlappingLinks() {

        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("contig0", 2000));
        dictionary.addSequence(new SAMSequenceRecord("contig1", 2000));
        final SVIntervalTree<EvidenceTargetLink> tree = new SVIntervalTree<>();

        final EvidenceTargetLink link = new EvidenceTargetLink(
                new StrandedInterval(new SVInterval(0, 100, 200), true),
                new StrandedInterval(new SVInterval(0, 600, 700), false),
                0, 0, Collections.emptySet(), Collections.emptySet());
        tree.put(new SVInterval(0, 100, 200), link);
        tree.put(new SVInterval(0, 600, 700), link);

        Assert.assertTrue(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 300, 400), tree, 50, dictionary).isEmpty());
        Assert.assertTrue(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 300, 400), tree, 100, dictionary).isEmpty());
        Assert.assertTrue(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(1, 710, 810), tree, 100, dictionary).isEmpty());
        Assert.assertTrue(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(1, 150, 250), tree, 0, dictionary).isEmpty());
        Assert.assertEquals(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 210, 310), tree, 100, dictionary), Collections.singletonList(link));
        Assert.assertEquals(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 300, 400), tree, 101, dictionary), Collections.singletonList(link));
        Assert.assertEquals(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 50, 150), tree, 0, dictionary), Collections.singletonList(link));
        Assert.assertEquals(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 550, 650), tree, 0, dictionary), Collections.singletonList(link));
        Assert.assertEquals(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 800, 900), tree, 50, dictionary), Collections.emptyList());
        Assert.assertEquals(LargeSimpleSVFactory.localOverlappingLinks(new SVInterval(0, 650, 750), tree, 0, dictionary), Collections.singletonList(link));
    }

    @DataProvider(name = "supportedByHMMData")
    public Object[][] getSupportedByHMMTestData() {
        return new Object[][]{
                {Arrays.asList(1.1, 1.02, 0.95, 1.03, 1.0, 0.87, 1.22, 1.05, 1.04, 1.37, 0.9, 1.15), false},
                {Arrays.asList(1.4, 1.5, 1.8, 1.6, 1.5, 1.4, 1.45, 1.65, 1.3, 1.6, 1.5), true},
                {Arrays.asList(0.1, 0.2, 0.1, 0.05, 0.01, 0.3, 0.2, 0.05, 0.01, 0.01, 0.25), false},
                {Arrays.asList(1.1, 0.9, 0.8, 0.9, 1.2, 1.7, 1.6, 1.4, 1.5, 1.5), true},
                {Arrays.asList(1.1, 0.9, 0.8, 0.9, 1.2, 1.1, 1.2, 1.0, 1.5, 1.5), false}
        };
    }

    @Test(groups = "sv",
            dataProvider = "supportedByHMMData")
    public void testSupportedByHMM(final List<Double> rawCopyRatios, final boolean expectedResult) {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("contig0", 2000));
        final List<CopyRatio> copyRatiosList = getCopyRatioBins(rawCopyRatios, "contig0", 0, 10);
        final LargeSimpleSVFactory factory = getEmptyTandemDuplicationFactory(Collections.emptyList(), Collections.emptyList(), Collections.emptyList(), dictionary);
        factory.arguments.hmmTransitionProb = 0.01;
        factory.arguments.hmmValidStatesMinFraction = 0.5;
        Assert.assertEquals(factory.supportedByHMM(copyRatiosList), expectedResult);
    }

    private List<CopyRatio> getCopyRatioBins(final List<Double> rawCopyRatios, final String contig, final int start, final int binSize) {
        final List<CopyRatio> copyRatiosList = new ArrayList(rawCopyRatios.size());
        for (int i = 0; i < rawCopyRatios.size(); i++) {
            copyRatiosList.add(new CopyRatio(new SimpleInterval(contig, start + 1 + i * binSize, start + (i + 1) * binSize), Math.log(rawCopyRatios.get(i)) / Math.log(2.0)));
        }
        return copyRatiosList;
    }
}