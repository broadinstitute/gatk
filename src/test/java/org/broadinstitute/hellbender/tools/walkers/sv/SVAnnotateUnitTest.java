package org.broadinstitute.hellbender.tools.walkers.sv;

import org.broadinstitute.hellbender.tools.spark.sv.utils.ClosedSVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

import static java.util.Objects.isNull;

public class SVAnnotateUnitTest {
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

    // transcription start sites and initTree() method and contigNames for testing annotateNearestTranscriptionStartSite()
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

    private final String[] contigNames = {"chr1", "chr2", "chr3", "chr4"};

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
}