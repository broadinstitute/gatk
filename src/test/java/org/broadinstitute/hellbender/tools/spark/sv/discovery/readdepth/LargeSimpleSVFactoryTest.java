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
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.StrandedInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;

public class LargeSimpleSVFactoryTest extends BaseTest {

    private static LargeSimpleSVFactory getEmptyTandemDuplicationFactory() {
        final SVIntervalTree<EvidenceTargetLink> intrachromosomalLinkTree = new SVIntervalTree<>();
        final SVIntervalTree<EvidenceTargetLink> interchromosomalLinkTree = new SVIntervalTree<>();
        final SVIntervalTree<GATKRead> contigTree = new SVIntervalTree<>();
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection arguments = new StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromReadDepthArgumentCollection();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();

        final SampleLocatableMetadata sampleLocatableMetadata = new SimpleSampleLocatableMetadata("test", dictionary);
        final CalledCopyRatioSegmentCollection calledCopyRatioSegmentCollection = new CalledCopyRatioSegmentCollection(sampleLocatableMetadata, Collections.emptyList());
        final CopyRatioCollection copyRatioSegmentCollection = new CopyRatioCollection(sampleLocatableMetadata, Collections.emptyList());

        final OverlapDetector<CalledCopyRatioSegment> copyRatioSegmentOverlapDetector = calledCopyRatioSegmentCollection.getOverlapDetector();
        final OverlapDetector<CopyRatio> copyRatioOverlapDetector = copyRatioSegmentCollection.getOverlapDetector();

        return new LargeTandemDuplicationFactory(intrachromosomalLinkTree, interchromosomalLinkTree, contigTree, arguments, copyRatioSegmentOverlapDetector, copyRatioOverlapDetector, dictionary);
    }

    @Test(groups = "sv")
    public void testCreate() {
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
}