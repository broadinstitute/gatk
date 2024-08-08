package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SVStratificationEngineUnitTest extends GATKBaseTest {

    private static final GATKPath CONFIG_FILE_PATH = new GATKPath(toolsTestDir + "/sv/sv_stratify_config.tsv");

    private static final String CONTEXT_1_NAME = "context1";
    private static final String CONTEXT_2_NAME = "context2";

    private static final List<Locatable> CONTEXT_1_INTERVALS = Lists.newArrayList(new SimpleInterval("chr1", 1000, 2000));
    private static final List<Locatable> CONTEXT_2_INTERVALS = Lists.newArrayList(new SimpleInterval("chr2", 1000, 2000));

    @Test
    public void testAddContext() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        engine.addContext(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        Assert.assertNotNull(engine.getContextIntervals(CONTEXT_1_NAME));
        Assert.assertNull(engine.getContextIntervals(CONTEXT_2_NAME));
    }
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddDuplicateContext() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        engine.addContext(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        engine.addContext(CONTEXT_1_NAME, CONTEXT_2_INTERVALS);
    }

    @Test
    public void testAddStratification() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        engine.addContext(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        engine.addStratification(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, CONTEXT_1_NAME);
        final Collection<SVStatificationEngine.SVStratification> stratificationCollection = engine.getStratifications();
        Assert.assertNotNull(stratificationCollection);
        Assert.assertEquals(stratificationCollection.size(), 1);
        final SVStatificationEngine.SVStratification stratification = stratificationCollection.iterator().next();
        Assert.assertNotNull(stratification);
        Assert.assertEquals(stratification.getSvType(), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertNotNull(stratification.getMinSize());
        Assert.assertEquals(stratification.getMinSize().intValue(), 50);
        Assert.assertNotNull(stratification.getMaxSize());
        Assert.assertEquals(stratification.getMaxSize().intValue(), 500);
        Assert.assertEquals(stratification.getContextName(), CONTEXT_1_NAME);
    }

    @Test
    public void testCreate() {
        final Map<String, List<Locatable>> map = new HashMap();
        map.put(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        map.put(CONTEXT_2_NAME, CONTEXT_2_INTERVALS);
        final SVStatificationEngine engine = SVStatificationEngine.create(map, CONFIG_FILE_PATH);
        Assert.assertNotNull(engine);
        Assert.assertNotNull(engine.getContextIntervals(CONTEXT_1_NAME));
        Assert.assertEquals(engine.getStratifications().size(), 9);
    }

    @DataProvider(name="testGetMatchVariantsData")
    public Object[][] testGetMatchVariantsData() {
        return new Object[][] {

                // DEL

                // Outside context interval
                { "chr1", 100, "chr1", 200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                { "chr1", 2000, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                // Simple match
                { "chr1", 1100, "chr1", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr1", 900, "chr1", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr1", 900, "chr1", 1900, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr1", 1100, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr1", 800, "chr1", 2100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr1", 999, "chr1", 2001, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr2", 1100, "chr2", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context2" },
                // Wrong contig
                { "chr3", 1100, "chr3", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                // Barely match
                { "chr1", 1000, "chr1", 3001, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr1", 2, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                { "chr1", 500, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                // Barely miss overlap threshold
                { "chr1", 1000, "chr1", 3002, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },
                // Barely large enough
                { "chr1", 1100, "chr1", 1149, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, "DEL_50_5000_context1" },
                // Too small
                { "chr1", 1100, "chr1", 1148, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, null, null },

                // INS

                // In context
                { "chr1", 1100, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 100, "INS_context1" },
                { "chr2", 1100, "chr2", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 300, "INS_300_400_context2" },
                { "chr2", 1100, "chr2", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 399, "INS_300_400_context2" },
                // SVLEN should not matter
                { "chr1", 1100, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 1, "INS_context1" },
                { "chr1", 1100, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 10000, "INS_context1" },
                // Out of context
                { "chr1", 100, "chr1", 100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 100, null },
                // Out of size range for context2
                { "chr2", 1100, "chr2", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 1000, null },
                { "chr2", 1100, "chr2", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS, 400, null },

                // BND

                // Both ends
                { "chr1", 1000, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, "BND_context1" },
                { "chr1", 2000, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, "BND_context1" },
                // One end only
                { "chr1", 500, "chr1", 900, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, null },
                { "chr1", 1500, "chr1", 3000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, null },
                // No ends
                { "chr1", 500, "chr1", 3000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND, null, null },
        };
    }

    @Test(dataProvider = "testGetMatchVariantsData")
    public void testGetMatchVariants(final String chromA, final int posA, final String chromB, final int posB,
                                     final GATKSVVCFConstants.StructuralVariantAnnotationType svType,
                                     final Integer svlen,
                                     final String expectedStratName) {
        final Map<String, List<Locatable>> map = new HashMap<>();
        map.put(CONTEXT_1_NAME, CONTEXT_1_INTERVALS);
        map.put(CONTEXT_2_NAME, CONTEXT_2_INTERVALS);
        final SVStatificationEngine engine = SVStatificationEngine.create(map, CONFIG_FILE_PATH);
        final SVCallRecord record;
        if (svType == GATKSVVCFConstants.StructuralVariantAnnotationType.INS) {
            record = SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates(chromA, posA, svlen);
        } else {
            record = SVTestUtils.newCallRecordWithCoordinatesAndType("record1", chromA, posA, chromB, posB, svType);
        }
        final SVStatificationEngine.SVStratification result = engine.getMatch(record, 0.5, 0, 2);
        if (expectedStratName == null) {
            Assert.assertNull(result);
        } else {
            Assert.assertNotNull(result);
            Assert.assertEquals(result.getName(), expectedStratName);
        }
    }

    @Test
    public void testTestAddStratificationInnerClass() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification stratification = engine.new SVStratification(GATKSVVCFConstants.StructuralVariantAnnotationType.DEL, 50, 500, CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS));
        engine.addStratification(stratification);
        final Collection<SVStatificationEngine.SVStratification> stratificationCollection = engine.getStratifications();
        Assert.assertNotNull(stratificationCollection);
        Assert.assertEquals(stratificationCollection.size(), 1);
        final SVStatificationEngine.SVStratification stratificationOut = stratificationCollection.iterator().next();
        Assert.assertNotNull(stratificationOut);
        Assert.assertEquals(stratificationOut.getSvType(), GATKSVVCFConstants.StructuralVariantAnnotationType.DEL);
        Assert.assertNotNull(stratificationOut.getMinSize());
        Assert.assertEquals(stratificationOut.getMinSize().intValue(), 50);
        Assert.assertNotNull(stratificationOut.getMaxSize());
        Assert.assertEquals(stratificationOut.getMaxSize().intValue(), 500);
        Assert.assertEquals(stratificationOut.getContextName(), CONTEXT_1_NAME);
    }

    @Test
    public void testStratificationMatches() {
    }

    @Test
    public void testMatchesType() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                100, 500,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );
        Assert.assertTrue(strat.matchesType(SVTestUtils.newCallRecordWithLengthAndType(null, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesType(SVTestUtils.newCallRecordWithLengthAndType(null, GATKSVVCFConstants.StructuralVariantAnnotationType.DUP)));
    }

    @Test
    public void testMatchesSizeSimple() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                100, 500,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(499, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(50, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(500, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeNoMin() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, 500,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(499, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(500, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeNoMax() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                50, null,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(49, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(Integer.MAX_VALUE - 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeNoMinOrMax() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.DEL,
                null, null,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(Integer.MAX_VALUE - 1, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL)));
    }

    @Test
    public void testMatchesSizeInsertion() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.INS,
                100, 500,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(100, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(499, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(50, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
        Assert.assertFalse(strat.matchesSize(SVTestUtils.newCallRecordWithLengthAndType(500, GATKSVVCFConstants.StructuralVariantAnnotationType.INS)));
    }

    @Test
    public void testMatchesSizeBnd() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, null,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );
        Assert.assertTrue(strat.matchesSize(SVTestUtils.newBndCallRecordWithStrands(true, false)));
    }

    @Test
    public void testMatchesContext() {
        final SVStatificationEngine engine = new SVStatificationEngine();
        final SVStatificationEngine.SVStratification strat = engine.new SVStratification(
                GATKSVVCFConstants.StructuralVariantAnnotationType.BND,
                null, null,
                CONTEXT_1_NAME, OverlapDetector.create(CONTEXT_1_INTERVALS)
        );

        // DEL
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1000, "chr1", 1500, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0.5, 0, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 500, "chr1", 1500, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0.5, 0, 1));
        Assert.assertFalse(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 499, "chr1", 1499, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0.5, 0, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 900, "chr1", 1300, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0.5, 0, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 900, "chr1", 1300, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0.5, 1, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1999, "chr1", 2000000, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0, 1, 1));
        Assert.assertFalse(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 500, "chr1", 600, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0, 2, 1));
        Assert.assertFalse(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 500, "chr1", 1100, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0, 2, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1100, "chr1", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                0, 2, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1100, "chr1", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.DEL),
                1, 2, 1));

        // INS
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates("chr1", 1100, 100),
                0.1, 0, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates("chr1", 1100, 100000),
                0.1, 0, 1));
        Assert.assertFalse(strat.matchesContext(SVTestUtils.newCallRecordInsertionWithLengthAndCoordinates("chr1", 999, 100),
                0.1, 0, 1));

        // BND
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1000, "chr1", 1200, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1000, "chr1", 3000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 500, "chr1", 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, 1));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1000, "chr1", 1999, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, 2));
        Assert.assertTrue(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1000, "chr1", 2000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, 2));
        Assert.assertFalse(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 1000, "chr1", 2001, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, 2));
        Assert.assertFalse(strat.matchesContext(SVTestUtils.newCallRecordWithCoordinatesAndType("record", "chr1", 999, "chr1", 1000, GATKSVVCFConstants.StructuralVariantAnnotationType.BND),
                0.5, 2, 2));
    }

    @Test
    public void testMatchesContextIntrachromosomal() {
    }

    @Test
    public void testMatchesContextOverlapFraction() {
    }

    @Test
    public void testMatchesContextBreakpointOverlap() {
    }

    @Test
    public void testCountAnyContextOverlap() {
    }

    @Test
    public void testIsMutuallyExclusive() {
    }

    @Test
    public void testGetSvType() {
    }

    @Test
    public void testGetMinSize() {
    }

    @Test
    public void testGetMaxSize() {
    }

    @Test
    public void testGetContextName() {
    }

    @Test
    public void testGetName() {
    }
}