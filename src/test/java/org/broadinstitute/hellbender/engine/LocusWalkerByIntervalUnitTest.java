package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class LocusWalkerByIntervalUnitTest extends CommandLineProgramTest {

    @CommandLineProgramProperties(
            summary = "Dummy that reads file and counts how many pileup elements are transformed",
            oneLineSummary = "none",
            programGroup = TestProgramGroup.class
    )
    private static class TestTransformedLocusWalker extends LocusWalkerByInterval {
        public Map<Locatable, Integer> overlapBases = new HashMap<>();
        public Map<Locatable, Boolean> wasOpened = new HashMap<>();
        public Map<Locatable, Boolean> wasClosed = new HashMap<>();

        public TestTransformedLocusWalker(Collection<Locatable> objectsToTestOverlapTo) {
            for (Locatable l : objectsToTestOverlapTo) {
                overlapBases.put(l, 0);
                wasOpened.put(l, false);
                wasClosed.put(l, false);
            }
        }

        @Override
        public List<Locatable> getIntervalObjectsToQueryOver() {
            return new ArrayList<>(overlapBases.keySet());
        }

        @Override
        public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext, Set<Locatable> activeIntervals) {
            for (Locatable l : activeIntervals) {
                overlapBases.put(l, overlapBases.get(l) + 1);
            }
        }

        @Override
        public Object onTraversalSuccess() {
            for (Locatable l : overlapBases.keySet()) {
                // Assert that every interval was closed that was opened
                Assert.assertEquals(wasOpened.get(l), wasClosed.get(l));
            }
            return null;
        }

        @Override
        public void onIntervalStart(Locatable activeInterval) { wasOpened.put(activeInterval, true);}
        @Override
        public void onIntervalEnd(Locatable activeInterval) { wasClosed.put(activeInterval, true);}
        public boolean emitEmptyLoci() {
            return true;
        }
    }

    @DataProvider
    public Object[][] getOverlapsAndOverlappingDataTestCases() {
        return new Object[][] {
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:1000-2000")},
                        new int[]{1001}},
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:500-2500")},
                        new int[]{1001}},
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("21:1000-2000")},
                        new int[]{0}},
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:1000-2000"), new SimpleInterval("20:1500-2000")},
                        new int[]{1001, 501}},
                {Arrays.asList("20:1000-2000", "20:2001-3000"),
                        new Locatable[]{new SimpleInterval("20:1000-3000"), new SimpleInterval("20:1500-2000")},
                        new int[]{2001, 501}},
                {Arrays.asList("20:1000-2000", "20:3001-4000"),
                        new Locatable[]{new SimpleInterval("20:1000-2500"), new SimpleInterval("20:1500-3500")},
                        new int[]{1001, 1001}},
                {Arrays.asList("20:1000-2000", "21:3001-4000"),
                        new Locatable[]{new SimpleInterval("20:1000-2500"), new SimpleInterval("20:100-500") , new SimpleInterval("21:100-500")},
                        new int[]{1001, 0, 0}},
                {Arrays.asList("20:1000-2000", "20:3001-4000"),
                        new Locatable[]{ArtificialReadUtils.createHeaderlessSamBackedRead("foo", "20", 1000, 500), VariantContextTestUtils.makeHomRef("20", 1300, 20, 1700), new SimpleInterval("20:1500-3500")},
                        new int[]{500, 401 , 1001}},


                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:1100-1250")},
                        new int[]{151}},
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:500-1200")},
                        new int[]{201}},
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:1550-2200")},
                        new int[]{451}},
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:400-500")},
                        new int[]{0}},
                {Arrays.asList("20:1000-2000"),
                        new Locatable[]{new SimpleInterval("20:2500-2600")},
                        new int[]{0}},
        };
    }


    @Test(dataProvider = "getOverlapsAndOverlappingDataTestCases")
    public void testOverlappingBasesCoverageInformation(List<String> inputIntervals, Locatable[] locatablesToQuery, int[] expectedApplyCounts) throws IOException {

        final TestTransformedLocusWalker tool = new TestTransformedLocusWalker(Arrays.asList(locatablesToQuery));
        String readInput = getTestDataDir() + "/../engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam";
        final ArrayList<String> args = new ArrayList<>();
        args.add("-I"); args.add(readInput);
        args.add("-R"); args.add(b37_reference_20_21);

        for (String interval : inputIntervals) {
            args.add("-L");
            args.add(interval);
        }

        tool.instanceMain(args.toArray(new String[args.size()]));

        // Assert that each of the overlapping was called in apply the expected number of itmes
        for (int i = 0; i < locatablesToQuery.length; i++) {
            Assert.assertEquals((int)tool.overlapBases.get(locatablesToQuery[i]), expectedApplyCounts[i]);
            if (expectedApplyCounts[i] > 0) {
                Assert.assertTrue(tool.wasOpened.get(locatablesToQuery[i]));
                Assert.assertTrue(tool.wasClosed.get(locatablesToQuery[i]));
            } else {
                Assert.assertFalse(tool.wasOpened.get(locatablesToQuery[i]));
                Assert.assertFalse(tool.wasClosed.get(locatablesToQuery[i]));
            }
        }

    }
}