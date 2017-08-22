package org.broadinstitute.hellbender.engine;

import htsjdk.variant.variantcontext.VariantContext;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class MultiVariantDataSourceUnitTest extends GATKBaseTest {
    private static final String ENGINE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final String MULTI_VARIANT_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/MultiVariantDataSource/";

    private static final File QUERY_TEST_VCF = new File(ENGINE_TEST_DIRECTORY + "feature_data_source_test.vcf");
    private static final File QUERY_TEST_GVCF = new File(ENGINE_TEST_DIRECTORY + "feature_data_source_test_gvcf.vcf");

    private static final File baseVariants = new File(MULTI_VARIANT_TEST_DIRECTORY, "baseVariants.vcf");
    private static final File baseVariantsAlternateDictionary = new File(MULTI_VARIANT_TEST_DIRECTORY, "baseVariantsAlternateDictionary.vcf");
    private static final File baseVariantsConflictingDictionary = new File(MULTI_VARIANT_TEST_DIRECTORY, "baseVariantsConflictingDictionary.vcf");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRequireFeatureInput() {
        new MultiVariantDataSource(Collections.emptyList(), FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testRejectNonExistentFile() {
        new MultiVariantDataSource(
                Collections.singletonList(
                        new FeatureInput<>(GATKBaseTest.getSafeNonExistentFile("nonexistent.vcf").getAbsolutePath(), "nonexistent")),
                FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRejectNullFile() {
        new MultiVariantDataSource(
                Collections.singletonList(
                        new FeatureInput<>(null, "sourceName1")),
                FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);
    }

    @Test(expectedExceptions = UserException.class)
    public void testQueryOverUnindexedFile() {
        try (MultiVariantDataSource multiVariantSource = new MultiVariantDataSource(
                Collections.singletonList(
                        new FeatureInput<>(new File(ENGINE_TEST_DIRECTORY + "unindexed.vcf").getAbsolutePath(), "unindexed")),
                FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            multiVariantSource.query(new SimpleInterval("1", 1, 1));
        }
    }

    @Test
    public void testGetName() {
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();
        featureInputs.add(new FeatureInput<>(baseVariants.getAbsolutePath(), "sourceName1"));

        try (final MultiVariantDataSource multiVariantSource =
                new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            Assert.assertTrue(multiVariantSource.getName().contains("sourceName1"));
        }

        featureInputs.add(new FeatureInput<>(baseVariantsAlternateDictionary.getAbsolutePath(), "sourceName2"));
        try (final MultiVariantDataSource multiVariantSource =
                     new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            String name = multiVariantSource.getName();
            Assert.assertTrue(name.contains("sourceName1"));
            Assert.assertTrue(name.contains("sourceName2"));
        }
    }

    @Test
    public void testGetSequenceDictionaryCompatible() {
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();

        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf").getAbsolutePath(), "interleavedVariants_1"));
        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf").getAbsolutePath(), "interleavedVariants_2"));

        try (final MultiVariantDataSource multiVariantSource =
                     new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            Assert.assertEquals(multiVariantSource.getSequenceDictionary().getSequences().size(), 4);
        }
    }

    @Test
    public void testGetSequenceDictionaryAlternate() {
        // tests the case where the files have alternate/disjoint contig sets
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();
        featureInputs.add(new FeatureInput<>(baseVariants.getAbsolutePath(), "baseVariants"));
        featureInputs.add(new FeatureInput<>(baseVariantsAlternateDictionary.getAbsolutePath(), "baseVariantsAlternateDictionary"));

        try (final MultiVariantDataSource multiVariantSource =
                     new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            Assert.assertEquals(multiVariantSource.getSequenceDictionary().getSequences().size(), 8);
        }
    }

    @Test(expectedExceptions = UserException.IncompatibleSequenceDictionaries.class)
    public void testGetSequenceDictionaryIncompatible() {
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();
        featureInputs.add(new FeatureInput<>(baseVariants.getAbsolutePath(), "baseVariants"));
        featureInputs.add(new FeatureInput<>(baseVariantsConflictingDictionary.getAbsolutePath(), "baseVariantsConflictingDictionary"));

        try (final MultiVariantDataSource multiVariantSource =
                     new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            multiVariantSource.getSequenceDictionary();
        }
    }

    @Test
    public void testIterator() {
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();

        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf").getAbsolutePath(), "interleavedVariants_1"));
        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf").getAbsolutePath(), "interleavedVariants_2"));

        try (MultiVariantDataSource multiVariantSource = new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            int count = 0;
            for (final VariantContext vc: multiVariantSource) {
                count++;
            };
            Assert.assertEquals(count, 26);
        }
    }

    @Test
    public void testIteratorOverlapping() {
        //Test interleaved files that include some variants that start at the same position in both files
        String expectedIDOrder[] = new String[] {
                "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n",
                "o", "o_overlap",
                "p", "q", "r", "s", "t", "u", "v", "w",
                "x", "x_overlap",
                "y", "z"
        };
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();

        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1_WithOverlap.vcf").getAbsolutePath(),
                "interleavedVariants_1_WithOverlap"));
        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2_WithOverlap.vcf").getAbsolutePath(),
                "interleavedVariants_2_WithOverlap"));

        try (final MultiVariantDataSource multiVariantSource =
                     new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            int count = 0;
            for (final VariantContext vc: multiVariantSource) {
                Assert.assertEquals(vc.getID(), expectedIDOrder[count]);
                count++;
            };
            Assert.assertEquals(count, 28);
        }
    }

    @Test
    public void testSerialQueries() {
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();

        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf").getAbsolutePath(), "interleavedVariants_1"));
        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf").getAbsolutePath(), "interleavedVariants_2"));

        try (final MultiVariantDataSource multiVariantSource =
                     new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            int count = 0;

            Iterator<VariantContext> it = multiVariantSource.query(new SimpleInterval("1", 1, 1200));
            while (it.hasNext()) {
                it.next();
                count++;
            };
            Assert.assertEquals(count, 14);

            count = 0;
            it = multiVariantSource.query(new SimpleInterval("2", 200, 600));
            while (it.hasNext()) {
                it.next();
                count++;
            };
            Assert.assertEquals(count, 3);
        }

    }

    @Test
    public void testSetIntervals() {
        List<FeatureInput<VariantContext>> featureInputs = new ArrayList<>();
        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_1.vcf").getAbsolutePath(), "interleavedVariants_1"));
        featureInputs.add(new FeatureInput<>(
                new File(MULTI_VARIANT_TEST_DIRECTORY, "interleavedVariants_2.vcf").getAbsolutePath(), "interleavedVariants_2"));

        try (final MultiVariantDataSource multiVariantSource =
                     new MultiVariantDataSource(featureInputs, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES)) {
            int count = 0;
            multiVariantSource.setIntervalsForTraversal(
                    Arrays.asList(
                            new SimpleInterval("1", 1, 1200),
                            new SimpleInterval("2", 200, 600)
                    )
            );
            for (final VariantContext vc: multiVariantSource) {
                count++;
            };
            Assert.assertEquals(count, 17);
        }
    }

    @DataProvider(name = "CompleteIterationTestData")
    public Object[][] getCompleteIterationTestData() {
        // File to iterate over + Expected Variant ID(s)
        return new Object[][] {
                { QUERY_TEST_VCF, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z") },
                // This file has no index and contains no embedded sequence dictionary, so its rejected.
                //{ UNINDEXED_VCF, Arrays.asList("a", "b", "c") }
        };
    }

    @Test(dataProvider = "CompleteIterationTestData")
    public void testCompleteIterationOverFile( final File vcfFile, final List<String> expectedVariantIDs ) {
        try ( MultiVariantDataSource multiVariantSource =
                  new MultiVariantDataSource(
                    Collections.singletonList(new FeatureInput<>(
                        vcfFile.getAbsolutePath(), vcfFile.getName())), FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES) ) {
            Iterator<VariantContext> iter = multiVariantSource.iterator();

            checkTraversalResults(iter, expectedVariantIDs, vcfFile, null);
        }
    }

    @DataProvider(name = "TraversalByIntervalsTestData")
    public Object[][] getTraversalByIntervalsTestData() {
        // Intervals for traversal + expected Variant IDs
        return new Object[][] {
                // Single interval
                { Arrays.asList(new SimpleInterval("1", 100, 200)), Arrays.asList("a", "b", "c") },

                // Two non-adjacent intervals on the same contig
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Some records overlap multiple intervals, and there are gaps between intervals
                { Arrays.asList(new SimpleInterval("1", 100, 203), new SimpleInterval("1", 205, 284), new SimpleInterval("1", 286, 1000)), Arrays.asList("a", "b", "c", "d", "e", "f", "h", "i", "j", "k") },

                // Some records overlap multiple intervals, and no gaps between intervals
                { Arrays.asList(new SimpleInterval("1", 100, 203), new SimpleInterval("1", 204, 285), new SimpleInterval("1", 286, 1000)), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k") },

                // Two intervals on different contigs
                { Arrays.asList(new SimpleInterval("1", 200, 300), new SimpleInterval("2", 500, 600)), Arrays.asList("b", "c", "d", "e", "f", "g", "h", "p", "q") },

                // More than two intervals spanning different contigs, and some records overlap multiple intervals
                { Arrays.asList(new SimpleInterval("1", 200, 203), new SimpleInterval("1", 205, 285), new SimpleInterval("2", 200, 548), new SimpleInterval("2", 550, 650), new SimpleInterval("4", 700, 800)), Arrays.asList("b", "c", "d", "e", "f", "g", "o", "p", "q", "r", "y", "z") },

                // One interval with no overlapping records at the beginning of interval list
                { Arrays.asList(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Multiple intervals with no overlapping records at the beginning of interval list
                { Arrays.asList(new SimpleInterval("1", 1, 50), new SimpleInterval("1", 60, 70), new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // One interval with no overlapping records in the middle of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 500, 600), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Multiple intervals with no overlapping records in the middle of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 500, 600), new SimpleInterval("1", 700, 800), new SimpleInterval("1", 1000, 2000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // One interval with no overlapping records at the end of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000), new SimpleInterval("1", 2000, 3000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // Multiple intervals with no overlapping records at the end of interval list
                { Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("1", 1000, 2000), new SimpleInterval("1", 2000, 3000), new SimpleInterval("1", 4000, 5000)), Arrays.asList("a", "b", "c", "j", "k", "l", "m", "n") },

                // No records overlap any intervals
                { Arrays.asList(new SimpleInterval("1", 1, 99), new SimpleInterval("1", 287, 290), new SimpleInterval("1", 500, 600), new SimpleInterval("2", 201, 524), new SimpleInterval("2", 1000, 2000), new SimpleInterval("4", 1, 500)), Collections.<String>emptyList() },

                // No intervals (should traverse the entire file)
                { Collections.<SimpleInterval>emptyList(), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z") }
        };
    }

    @Test(dataProvider = "TraversalByIntervalsTestData")
    public void testTraversalByIntervals( final List<SimpleInterval> intervalsForTraversal, final List<String> expectedVariantIDs ) {
        try ( MultiVariantDataSource multiVariantSource =
              new MultiVariantDataSource(
                      Collections.singletonList(new FeatureInput<>(QUERY_TEST_VCF.getAbsolutePath(), QUERY_TEST_VCF.getName())),
                      FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES) ) {
            multiVariantSource.setIntervalsForTraversal(intervalsForTraversal);
            Iterator<VariantContext> iter = multiVariantSource.iterator();

            checkTraversalResults(iter, expectedVariantIDs, QUERY_TEST_VCF, intervalsForTraversal);
        }
    }

    private void checkTraversalResults( final Iterator<VariantContext> traversalResults, final List<String> expectedVariantIDs, final File vcfFile, final List<SimpleInterval> traversalIntervals ) {
        final String intervalString = traversalIntervals != null ? " with intervals " + traversalIntervals : "";

        int recordCount = 0;
        while ( traversalResults.hasNext() ) {
            VariantContext record = traversalResults.next();
            Assert.assertTrue(recordCount < expectedVariantIDs.size(), "Too many records returned during iteration over " + vcfFile.getAbsolutePath() + intervalString);
            Assert.assertEquals(record.getID(), expectedVariantIDs.get(recordCount),
                    "Record #" + (recordCount + 1) + " encountered in iteration over " + vcfFile.getAbsolutePath() + intervalString + " is incorrect");
            ++recordCount;
        }

        Assert.assertEquals(recordCount, expectedVariantIDs.size(), "Wrong number of records returned in iteration over " + vcfFile.getAbsolutePath() + intervalString);
    }

    @DataProvider(name = "GVCFQueryTestData")
    public Object[][] getGVCFQueryTestData() {

        // Ensure that queries on a FeatureDataSource take GVCF blocks into account when computing overlap.

        // Query interval + expected variant ID(s)
        return new Object[][] {
                { new SimpleInterval("1", 1, 99), Collections.<String>emptyList() },
                { new SimpleInterval("1", 50, 100), Arrays.asList("aa") },
                { new SimpleInterval("1", 50, 150), Arrays.asList("aa") },
                { new SimpleInterval("1", 100, 100), Arrays.asList("aa") },
                { new SimpleInterval("1", 100, 150), Arrays.asList("aa") },
                { new SimpleInterval("1", 100, 200), Arrays.asList("aa") },
                { new SimpleInterval("1", 150, 200), Arrays.asList("aa") },
                { new SimpleInterval("1", 150, 250), Arrays.asList("aa") },
                { new SimpleInterval("1", 200, 201), Arrays.asList("aa") },
                { new SimpleInterval("1", 201, 3000), Collections.<String>emptyList() }
        };
    }

    @Test(dataProvider = "GVCFQueryTestData")
    public void testQueryGVCF( final SimpleInterval queryInterval, final List<String> expectedVariantIDs ) {
        try ( MultiVariantDataSource multiVariantSource =
               new MultiVariantDataSource(
                       Collections.singletonList(new FeatureInput<>(QUERY_TEST_GVCF.getAbsolutePath(), QUERY_TEST_GVCF.getName())),
                       FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES) ) {
            Iterator<VariantContext> it = multiVariantSource.query(queryInterval);
            checkTraversalResults(it, expectedVariantIDs, QUERY_TEST_GVCF, Collections.singletonList(queryInterval));
        }
    }

}


