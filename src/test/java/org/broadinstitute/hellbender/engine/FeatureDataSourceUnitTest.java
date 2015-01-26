package org.broadinstitute.hellbender.engine;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class FeatureDataSourceUnitTest extends BaseTest {
    private static final String FEATURE_DATA_SOURCE_TEST_DIRECTORY = publicTestDir + "org/broadinstitute/hellbender/engine/";
    private static final File QUERY_TEST_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "feature_data_source_test.vcf");
    private static final File UNINDEXED_VCF = new File(FEATURE_DATA_SOURCE_TEST_DIRECTORY + "unindexed.vcf");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullFile() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(null, new VCFCodec());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleNullCodec() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, null);
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testHandleNonExistentFile() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(new File("/foo/bar/nonexistent.vcf"), new VCFCodec());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testHandleInvalidQueryLookahead() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec(), "MyName", -1);
    }

    @Test(expectedExceptions = UserException.class)
    public void testHandleQueryOverUnindexedFile() {
        try ( FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(UNINDEXED_VCF, new VCFCodec()) ) {
            featureSource.query(hg19GenomeLocParser.createGenomeLoc("1", 1, 1));  // Should throw, since we have no index
        }
    }

    @Test
    public void testGetCodecClass() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec());
        Assert.assertEquals(featureSource.getCodecClass(), VCFCodec.class, "Wrong codec class returned from getCodecClass()");
        featureSource.close();
    }

    @Test
    public void testGetFeatureType() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec());
        Assert.assertEquals(featureSource.getFeatureType(), VariantContext.class, "Wrong feature type returned from getFeatureType()");
        featureSource.close();
    }

    @Test
    public void testGetName() {
        FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec(), "CustomName");
        Assert.assertEquals(featureSource.getName(), "CustomName", "Wrong name returned from getName()");
        featureSource.close();
    }

    @DataProvider(name = "CompleteIterationTestData")
    public Object[][] getCompleteIterationTestData() {
        // File to iterate over + Expected Variant ID(s)
        return new Object[][] {
                { QUERY_TEST_VCF, Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z") },
                { UNINDEXED_VCF, Arrays.asList("a", "b", "c") }
        };
    }

    @Test(dataProvider = "CompleteIterationTestData")
    public void testCompleteIterationOverFile( final File vcfFile, final List<String> expectedVariantIDs ) {
        try ( FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(vcfFile, new VCFCodec()) ) {
            Iterator<VariantContext> iter = featureSource.iterator();

            int recordCount = 0;
            while ( iter.hasNext() ) {
                VariantContext record = iter.next();
                Assert.assertTrue(recordCount < expectedVariantIDs.size(), "Too many records returned during complete iteration over " + vcfFile.getAbsolutePath());
                Assert.assertEquals(record.getID(), expectedVariantIDs.get(recordCount),
                                    "Record #" + (recordCount + 1) + " encountered in iteration over " + vcfFile.getAbsolutePath() + " is incorrect");
                ++recordCount;
            }

            Assert.assertEquals(recordCount, expectedVariantIDs.size(), "Wrong number of records returned in complete iteration over " + vcfFile.getAbsolutePath());
        }
    }

    @DataProvider(name = "IndependentFeatureQueryTestData")
    public Object[][] getIndependentFeatureQueryTestData() {
        final GenomeLocParser parser = hg19GenomeLocParser;

        // Query Interval + Expected Variant ID(s)
        return new Object[][] {
                { parser.createGenomeLoc("1", 1, 99), Collections.<String>emptyList() },
                { parser.createGenomeLoc("1", 100, 100), Arrays.asList("a") },
                { parser.createGenomeLoc("1", 100, 200), Arrays.asList("a", "b", "c") },
                { parser.createGenomeLoc("1", 200, 202), Arrays.asList("b", "c") },
                { parser.createGenomeLoc("1", 200, 203), Arrays.asList("b", "c", "d") },
                { parser.createGenomeLoc("1", 201, 203), Arrays.asList("d") },
                { parser.createGenomeLoc("1", 204, 204), Arrays.asList("d") },
                { parser.createGenomeLoc("1", 206, 206), Arrays.asList("d") },
                { parser.createGenomeLoc("1", 207, 207), Collections.<String>emptyList() },
                { parser.createGenomeLoc("1", 200, 300), Arrays.asList("b", "c", "d", "e", "f", "g", "h") },
                { parser.createGenomeLoc("1", 275, 300), Arrays.asList("e", "f", "g", "h") },
                { parser.createGenomeLoc("1", 275, 284), Arrays.asList("e", "f") },
                { parser.createGenomeLoc("1", 284, 284), Arrays.asList("f") },
                { parser.createGenomeLoc("1", 284, 285), Arrays.asList("f", "g") },
                { parser.createGenomeLoc("1", 284, 286), Arrays.asList("f", "g", "h") },
                { parser.createGenomeLoc("1", 286, 286), Arrays.asList("f", "h") },
                { parser.createGenomeLoc("1", 287, 290), Collections.<String>emptyList() },
                { parser.createGenomeLoc("1", 999, 1000), Arrays.asList("i", "j", "k") },
                { parser.createGenomeLoc("1", 1000, 1001), Arrays.asList("j", "k") },
                { parser.createGenomeLoc("1", 1002, 1005), Arrays.asList("k") },
                { parser.createGenomeLoc("1", 1005, 1010), Collections.<String>emptyList() },
                { parser.createGenomeLoc("1", 1075, 1175), Arrays.asList("l", "m") },
                { parser.createGenomeLoc("1", 1075, 1176), Arrays.asList("l", "m", "n") },
                { parser.createGenomeLoc("1", 1077, 1176), Arrays.asList("m", "n") },
                { parser.createGenomeLoc("1", 1170, 1180), Arrays.asList("n") },
                { parser.createGenomeLoc("1", 1003, 1175), Arrays.asList("k", "l", "m") },
                { parser.createGenomeLoc("1", 1, 2000), Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n") },
                { parser.createGenomeLoc("1", 2000, 3000), Collections.<String>emptyList() },
                { parser.createGenomeLoc("1", 286, 1175), Arrays.asList("f", "h", "i", "j", "k", "l", "m") },
                { parser.createGenomeLoc("2", 200, 700), Arrays.asList("o", "p", "q", "r", "s") },
                { parser.createGenomeLoc("2", 201, 699), Arrays.asList("p", "q", "r") },
                { parser.createGenomeLoc("2", 550, 560), Arrays.asList("q") },
                { parser.createGenomeLoc("2", 549, 699), Arrays.asList("q", "r") },
                { parser.createGenomeLoc("2", 600, 700), Arrays.asList("r", "s") },
                { parser.createGenomeLoc("2", 701, 800), Collections.<String>emptyList() },
                { parser.createGenomeLoc("3", 1, 300), Arrays.asList("t", "u", "v") },
                { parser.createGenomeLoc("3", 300, 400), Arrays.asList("u", "v", "w") },
                { parser.createGenomeLoc("3", 301, 400), Arrays.asList("v", "w") },
                { parser.createGenomeLoc("3", 301, 399), Arrays.asList("v") },
                { parser.createGenomeLoc("3", 305, 400), Arrays.asList("w") },
                { parser.createGenomeLoc("3", 500, 600), Collections.<String>emptyList() },
                { parser.createGenomeLoc("4", 1, 1000), Arrays.asList("x", "y", "z") },
                { parser.createGenomeLoc("4", 600, 775), Arrays.asList("x", "y") },
                { parser.createGenomeLoc("4", 775, 776), Arrays.asList("y", "z") },
                { parser.createGenomeLoc("4", 777, 780), Arrays.asList("z") }
        };
    }

    /**
     * Tests correctness of the raw query operation (without any caching coming into play) by creating
     * a new FeatureDataSource for each query to eliminate caching as a factor.
     *
     * (We test the use case of re-using a single FeatureDataSource across queries in separate tests
     *  below to make sure the caching is working as expected)
     */
    @Test(dataProvider = "IndependentFeatureQueryTestData")
    public void testIndependentFeatureQuerying( final GenomeLoc queryInterval, final List<String> expectedVariantIDs ) {
        final FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec());

        // Use query() here rather than queryAndPrefetch() so that query() will have test coverage
        // (the other tests below use queryAndPrefetch())
        Iterator<VariantContext> featureIterator = featureSource.query(queryInterval);
        List<VariantContext> queryResults = new ArrayList<>();
        while ( featureIterator.hasNext() ) {
            queryResults.add(featureIterator.next());
        }
        featureSource.close();

        checkVariantQueryResults(queryResults, expectedVariantIDs, queryInterval);
    }

    private void checkVariantQueryResults( final List<VariantContext> queryResults, final List<String> expectedVariantIDs, final GenomeLoc queryInterval ) {
        Assert.assertEquals(queryResults.size(), expectedVariantIDs.size(), "Wrong number of records returned for query on interval " + queryInterval);
        for ( int recordIndex = 0; recordIndex < queryResults.size(); ++recordIndex ) {
            Assert.assertEquals(queryResults.get(recordIndex).getID(), expectedVariantIDs.get(recordIndex),
                                "Record #" + (recordIndex + 1) + " returned from query on interval " + queryInterval + " is incorrect");
        }
    }

    /**
     * This data provider serves up sets of multiple queries paired with a single FeatureDataSource.
     * The intention is to test caching behavior, including recovery from cache misses, and correctness
     * of the code that subsets the Feature cache down to Features overlapping the desired interval
     * (whereas the primary intent of the independent querying tests above was to test correctness of the
     * raw query operation without any caching involved).
     */
    @DataProvider(name = "SingleDataSourceMultipleQueriesTestData")
    @SuppressWarnings("unchecked")
    public Object[][] getSingleDataSourceMultipleQueriesTestData() {
        final GenomeLocParser parser = hg19GenomeLocParser;

        // Query set #1:
        // Re-use the queries + expected results from the IndependentFeatureQueryTestData DataProvider above,
        // but this time aggregated together and executed on the same FeatureDataSource
        List<Pair<GenomeLoc, List<String>>> aggregatedIndependentQueries = new ArrayList<>();
        Object[][] independentQueryTestData = getIndependentFeatureQueryTestData();
        for ( Object[] queryTest : independentQueryTestData ) {
            aggregatedIndependentQueries.add(Pair.of((GenomeLoc)queryTest[0], (List<String>)queryTest[1]));
        }

        // Query set #2:
        // Large intervals with regularly-increasing start positions. Represents typical query access patterns
        // ideal for the caching implementation in FeatureDataSource, as it minimizes cache misses.
        List<Pair<GenomeLoc, List<String>>> regularlyIncreasingQueries = Arrays.asList(
                Pair.of(parser.createGenomeLoc("1", 1, 100), Arrays.asList("a")),
                Pair.of(parser.createGenomeLoc("1", 50, 150), Arrays.asList("a")),
                Pair.of(parser.createGenomeLoc("1", 100, 200), Arrays.asList("a", "b", "c")),
                Pair.of(parser.createGenomeLoc("1", 150, 250), Arrays.asList("b", "c", "d")),
                Pair.of(parser.createGenomeLoc("1", 200, 300), Arrays.asList("b", "c", "d", "e", "f", "g", "h")),
                Pair.of(parser.createGenomeLoc("1", 250, 350), Arrays.asList("e", "f", "g", "h")),
                Pair.of(parser.createGenomeLoc("1", 300, 400), Collections.<String>emptyList()),
                Pair.of(parser.createGenomeLoc("1", 350, 450), Collections.<String>emptyList()),
                Pair.of(parser.createGenomeLoc("1", 950, 1050), Arrays.asList("i", "j", "k")),
                Pair.of(parser.createGenomeLoc("1", 1000, 1100), Arrays.asList("j", "k", "l")),
                // First cache miss here given the default queryLookaheadBases value of 1000
                Pair.of(parser.createGenomeLoc("1", 1050, 1150), Arrays.asList("l", "m")),
                Pair.of(parser.createGenomeLoc("1", 1100, 1200), Arrays.asList("m", "n")),
                Pair.of(parser.createGenomeLoc("1", 1150, 1250), Arrays.asList("m", "n")),
                Pair.of(parser.createGenomeLoc("1" ,1200, 1300), Collections.<String>emptyList()),
                // Second cache miss here as we change contigs
                Pair.of(parser.createGenomeLoc("2", 1, 100), Collections.<String>emptyList()),
                Pair.of(parser.createGenomeLoc("2", 100, 200), Arrays.asList("o")),
                Pair.of(parser.createGenomeLoc("2", 500, 600), Arrays.asList("p", "q")),
                Pair.of(parser.createGenomeLoc("2", 550, 650), Arrays.asList("q", "r")),
                Pair.of(parser.createGenomeLoc("2", 600, 700), Arrays.asList("r", "s")),
                Pair.of(parser.createGenomeLoc("2", 650, 750), Arrays.asList("s")),
                // Third cache miss (changing contigs again)
                Pair.of(parser.createGenomeLoc("3", 1, 200), Arrays.asList("t")),
                Pair.of(parser.createGenomeLoc("3", 300, 400), Arrays.asList("u", "v", "w")),
                Pair.of(parser.createGenomeLoc("3", 302, 350), Arrays.asList("v"))
        );

        // Query set #3:
        // Cache miss hell: lots of cache misses in this query set due to either backing up or skipping ahead too far
        List<Pair<GenomeLoc, List<String>>> cacheMissQueries = Arrays.asList(
                Pair.of(parser.createGenomeLoc("1", 100, 200), Arrays.asList("a", "b", "c")),
                // Cache miss due to backup
                Pair.of(parser.createGenomeLoc("1", 99, 205), Arrays.asList("a", "b", "c", "d")),
                // Cache miss due to jumping past query lookahead
                Pair.of(parser.createGenomeLoc("1", 2000, 3000), Collections.<String>emptyList()),
                // Cache miss due to backup
                Pair.of(parser.createGenomeLoc("1", 205, 285), Arrays.asList("d", "e", "f", "g")),
                Pair.of(parser.createGenomeLoc("1", 286, 400), Arrays.asList("f", "h")),
                // Cache miss due to contig change
                Pair.of(parser.createGenomeLoc("2", 200, 600), Arrays.asList("o", "p", "q")),
                // Cache miss due to backup
                Pair.of(parser.createGenomeLoc("2", 1, 100), Collections.<String>emptyList()),
                // Cache miss due to jumping past query lookahead
                Pair.of(parser.createGenomeLoc("2", 3000, 4000), Collections.<String>emptyList()),
                // Cache miss due to contig change
                Pair.of(parser.createGenomeLoc("1", 200, 300), Arrays.asList("b", "c", "d", "e", "f", "g", "h")),
                // Cache miss due to backup
                Pair.of(parser.createGenomeLoc("1", 100, 200), Arrays.asList("a", "b", "c")),
                // Cache miss due to backup
                Pair.of(parser.createGenomeLoc("1", 1, 1), Collections.<String>emptyList()),
                // Cache miss due to jumping past query lookahead
                Pair.of(parser.createGenomeLoc("1", 1100, 1200), Arrays.asList("m", "n"))
        );

        return new Object[][] {
                { aggregatedIndependentQueries },
                { regularlyIncreasingQueries },
                { cacheMissQueries }
        };
    }

    /**
     * Tests correctness of Feature caching behavior by executing multiple queries on the same FeatureDataSource
     */
    @Test(dataProvider = "SingleDataSourceMultipleQueriesTestData")
    public void testSingleDataSourceMultipleQueries( final List<Pair<GenomeLoc, List<String>>> testQueries ) {
        final FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(QUERY_TEST_VCF, new VCFCodec());

        // This test re-uses the same FeatureDataSource across queries to test caching of query results.
        for ( Pair<GenomeLoc, List<String>> testQuery : testQueries ) {
            final GenomeLoc queryInterval = testQuery.getLeft();
            final List<String> expectedVariantIDs = testQuery.getRight();

            final List<VariantContext> queryResults = featureSource.queryAndPrefetch(queryInterval);
            checkVariantQueryResults(queryResults, expectedVariantIDs, queryInterval);
        }

        featureSource.close();
    }

    /**************************************************
     * Direct testing on the FeatureCache inner class
     **************************************************/

    @SuppressWarnings("overrides")  // because I don't want to implement hashCode() but do need an equals() here
    private static class ArtificialTestFeature implements Feature {
        private String chr;
        private int start;
        private int end;

        public ArtificialTestFeature( final String chr, final int start, final int end ) {
            this.chr = chr;
            this.start = start;
            this.end = end;
        }

        @Override
        public String getChr() {
            return chr;
        }

        @Override
        public int getStart() {
            return start;
        }

        @Override
        public int getEnd() {
            return end;
        }

        @Override
        public boolean equals( Object other ) {
            if ( other == null || ! (other instanceof ArtificialTestFeature) ) {
                return false;
            }

            ArtificialTestFeature otherFeature = (ArtificialTestFeature)other;
            return chr.equals(otherFeature.getChr()) && start == otherFeature.getStart() && end == otherFeature.getEnd();
        }

        @Override
        public String toString() {
            return chr + ":" + start + "-" + end;   // (to improve output on test failures involving this class)
        }
    }

    private FeatureDataSource.FeatureCache<ArtificialTestFeature> initializeFeatureCache( final List<ArtificialTestFeature> features, final String cacheContig, final int cacheStart, final int cacheStop ) {
        GenomeLocParser parser = hg19GenomeLocParser;
        FeatureDataSource.FeatureCache<ArtificialTestFeature> cache = new FeatureDataSource.FeatureCache<>();

        cache.fill(features.iterator(), parser.createGenomeLoc(cacheContig, cacheStart, cacheStop));
        return cache;
    }

    @DataProvider(name = "FeatureCacheFillDataProvider")
    public Object[][] getFeatureCacheFillData() {
        return new Object[][] {
                { Arrays.asList(new ArtificialTestFeature("1", 1, 100), new ArtificialTestFeature("1", 50, 150),
                               new ArtificialTestFeature("1", 200, 300), new ArtificialTestFeature("1", 350, 400)),
                  "1", 1, 400 },
                { Arrays.asList(new ArtificialTestFeature("1", 1, 100)), "1", 1, 100 },
                { Collections.<ArtificialTestFeature>emptyList(), "1", 1, 1 }
        };
    }

    @Test(dataProvider = "FeatureCacheFillDataProvider")
    public void testCacheFill( final List<ArtificialTestFeature> features, final String cacheContig, final int cacheStart, final int cacheStop) {
        FeatureDataSource.FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(features, cacheContig, cacheStart, cacheStop);

        List<ArtificialTestFeature> cachedFeatures = cache.getCachedFeaturesUpToStopPosition(cacheStop);
        Assert.assertEquals(cache.getContig(), cacheContig, "Wrong contig reported by cache after fill");
        Assert.assertEquals(cache.getCacheStart(), cacheStart, "Wrong start position reported by cache after fill");
        Assert.assertEquals(cache.getCacheStop(), cacheStop, "Wrong stop position reported by cache after fill");
        Assert.assertEquals(cachedFeatures, features, "Wrong Features in cache after fill()");
    }

    @DataProvider(name = "FeatureCacheHitDetectionDataProvider")
    public Object[][] getFeatureCacheHitDetectionData() {
        List<ArtificialTestFeature> features = Arrays.asList(new ArtificialTestFeature("1", 1, 100),
                                                             new ArtificialTestFeature("1", 50, 150),
                                                             new ArtificialTestFeature("1", 200, 300));
        FeatureDataSource.FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(features, "1", 50, 250);
        GenomeLocParser parser = hg19GenomeLocParser;

        return new Object[][] {
                // Exact match for cache boundaries
                { cache, parser.createGenomeLoc("1", 50, 250), true },
                // Interval completely contained within cache boundaries
                { cache, parser.createGenomeLoc("1", 100, 200), true },
                // Interval left-aligned with cache boundaries
                { cache, parser.createGenomeLoc("1", 50, 100), true },
                // Interval right-aligned with cache boundaries
                { cache, parser.createGenomeLoc("1", 200, 250), true },
                // Interval overlaps, but is off the left edge of cache boundaries
                { cache, parser.createGenomeLoc("1", 49, 100), false },
                // Interval overlaps, but is off the right edge of cache boundaries
                { cache, parser.createGenomeLoc("1", 200, 251), false },
                // Interval does not overlap and is to the left of cache boundaries
                { cache, parser.createGenomeLoc("1", 1, 40), false },
                // Interval does not overlap and is to the right of cache boundaries
                { cache, parser.createGenomeLoc("1", 300, 350), false },
                // Interval is on different contig
                { cache, parser.createGenomeLoc("2", 50, 250), false }
        };
    }

    @Test(dataProvider = "FeatureCacheHitDetectionDataProvider")
    public void testCacheHitDetection( final FeatureDataSource.FeatureCache<ArtificialTestFeature> cache,
                                       final GenomeLoc testInterval, final boolean cacheHitExpectedResult ) {
        Assert.assertEquals(cache.cacheHit(testInterval), cacheHitExpectedResult,
                            "Cache hit detection failed for interval " + testInterval);
    }

    @DataProvider(name = "FeatureCacheTrimmingDataProvider")
    public Object[][] getFeatureCacheTrimmingData() {
        // Features are required to always be sorted by start position, but stop positions need not be sorted.
        // This complicates cache trimming.
        List<ArtificialTestFeature> feats = Arrays.asList(
                new ArtificialTestFeature("1", 1, 1),     // Feature 0
                new ArtificialTestFeature("1", 1, 100),   // Feature 1
                new ArtificialTestFeature("1", 1, 1),     // Feature 2
                new ArtificialTestFeature("1", 1, 50),    // Feature 3
                new ArtificialTestFeature("1", 1, 3),     // Feature 4
                new ArtificialTestFeature("1", 1, 5),     // Feature 5
                new ArtificialTestFeature("1", 5, 5),     // Feature 6
                new ArtificialTestFeature("1", 5, 50),    // Feature 7
                new ArtificialTestFeature("1", 5, 10),    // Feature 8
                new ArtificialTestFeature("1", 50, 100),  // Feature 9
                new ArtificialTestFeature("1", 50, 50),   // Feature 10
                new ArtificialTestFeature("1", 50, 200),  // Feature 11
                new ArtificialTestFeature("1", 100, 100), // Feature 12
                new ArtificialTestFeature("1", 100, 110), // Feature 13
                new ArtificialTestFeature("1", 100, 200), // Feature 14
                new ArtificialTestFeature("1", 100, 150), // Feature 15
                new ArtificialTestFeature("1", 100, 199)  // Feature 16
        );
        FeatureDataSource.FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(feats, "1", 1, 200);

        // Pairing of start position to which to trim the cache with the List of Features we expect to see
        // in the cache after trimming
        List<Pair<Long, List<ArtificialTestFeature>>> trimOperations = Arrays.asList(
                Pair.of(1l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(2l, Arrays.asList(feats.get(1), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(3l, Arrays.asList(feats.get(1), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(4l, Arrays.asList(feats.get(1), feats.get(3), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(5l, Arrays.asList(feats.get(1), feats.get(3), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(6l, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(10l, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(11l, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(50l, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(51l, Arrays.asList(feats.get(1), feats.get(9), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(100l, Arrays.asList(feats.get(1), feats.get(9), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(101l, Arrays.asList(feats.get(11), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(111l, Arrays.asList(feats.get(11), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(151l, Arrays.asList(feats.get(11), feats.get(14), feats.get(16))),
                Pair.of(151l, Arrays.asList(feats.get(11), feats.get(14), feats.get(16))),
                Pair.of(200l, Arrays.asList(feats.get(11), feats.get(14)))
        );

        return new Object[][] {
                { cache, trimOperations }
        };
    }

    @Test(dataProvider = "FeatureCacheTrimmingDataProvider")
    public void testCacheTrimming( final FeatureDataSource.FeatureCache<ArtificialTestFeature> cache, final List<Pair<Long, List<ArtificialTestFeature>>> trimOperations ) {
        // Repeatedly trim the cache to ever-increasing start positions, and verify after each trim operation
        // that the cache holds the correct Features in the correc order
        for ( Pair<Long, List<ArtificialTestFeature>> trimOperation : trimOperations ) {
            final long trimPosition = trimOperation.getLeft();
            final List<ArtificialTestFeature> expectedFeatures = trimOperation.getRight();

            cache.trimToNewStartPosition(trimPosition);

            final List<ArtificialTestFeature> actualFeatures = cache.getCachedFeaturesUpToStopPosition(cache.getCacheStop());
            Assert.assertEquals(actualFeatures, expectedFeatures, "Wrong Features in cache after trimming start position to " + trimPosition);
        }
    }

    @DataProvider(name = "FeatureCacheRetrievalDataProvider")
    public Object[][] getFeatureCacheRetrievalData() {
        List<ArtificialTestFeature> feats = Arrays.asList(
                new ArtificialTestFeature("1", 1, 1),      // Feature 0
                new ArtificialTestFeature("1", 1, 100),    // Feature 1
                new ArtificialTestFeature("1", 5, 5),      // Feature 2
                new ArtificialTestFeature("1", 10, 10),    // Feature 3
                new ArtificialTestFeature("1", 10, 100),   // Feature 4
                new ArtificialTestFeature("1", 50, 50),    // Feature 5
                new ArtificialTestFeature("1", 51, 55),    // Feature 6
                new ArtificialTestFeature("1", 52, 52),    // Feature 7
                new ArtificialTestFeature("1", 55, 60),    // Feature 8
                new ArtificialTestFeature("1", 75, 75),    // Feature 9
                new ArtificialTestFeature("1", 80, 100)    // Feature 10
        );
        FeatureDataSource.FeatureCache<ArtificialTestFeature> cache = initializeFeatureCache(feats, "1", 1, 100);

        // Pairing of end position with which to bound cache retrieval with the List of Features we expect to see
        // after retrieval
        List<Pair<Long, List<ArtificialTestFeature>>> retrievalOperations = Arrays.asList(
                Pair.of(100l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(80l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(79l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9))),
                Pair.of(80l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(75l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9))),
                Pair.of(74l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8))),
                Pair.of(54l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7))),
                Pair.of(52l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7))),
                Pair.of(51l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6))),
                Pair.of(50l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5))),
                Pair.of(49l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4))),
                Pair.of(10l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4))),
                Pair.of(9l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2))),
                Pair.of(5l, Arrays.asList(feats.get(0), feats.get(1), feats.get(2))),
                Pair.of(4l, Arrays.asList(feats.get(0), feats.get(1))),
                Pair.of(1l, Arrays.asList(feats.get(0), feats.get(1)))
        );

        return new Object[][] {
                { cache, retrievalOperations }
        };
    }

    @Test(dataProvider = "FeatureCacheRetrievalDataProvider")
    public void testCacheFeatureRetrieval( final FeatureDataSource.FeatureCache<ArtificialTestFeature> cache, final List<Pair<Long, List<ArtificialTestFeature>>> retrievalOperations ) {
        for ( Pair<Long, List<ArtificialTestFeature>> retrievalOperation: retrievalOperations ) {
            final long stopPosition = retrievalOperation.getLeft();
            final List<ArtificialTestFeature> expectedFeatures = retrievalOperation.getRight();

            final List<ArtificialTestFeature> actualFeatures = cache.getCachedFeaturesUpToStopPosition(stopPosition);
            Assert.assertEquals(actualFeatures, expectedFeatures, "Wrong Features returned in retrieval operation with stop position " + stopPosition);
        }
    }

    /**
     * Test caching a region with no Features. This should work (we should avoid going to disk
     * to look for new records when querying within such a region).
     */
    @Test
    public void testHandleCachingOfEmptyRegion() {
        GenomeLocParser parser = hg19GenomeLocParser;
        FeatureDataSource.FeatureCache<ArtificialTestFeature> cache = new FeatureDataSource.FeatureCache<>();
        List<ArtificialTestFeature> emptyRegion = new ArrayList<>();

        cache.fill(emptyRegion.iterator(), parser.createGenomeLoc("1", 1, 100));

        Assert.assertTrue(cache.isEmpty(), "Cache should be empty");
        Assert.assertTrue(cache.cacheHit(parser.createGenomeLoc("1", 1, 100)), "Unexpected cache miss");
        Assert.assertTrue(cache.cacheHit(parser.createGenomeLoc("1", 2, 99)), "Unexpected cache miss");

        Assert.assertEquals(cache.getCachedFeaturesUpToStopPosition(100), emptyRegion, "Should get back empty List for empty region");
        cache.trimToNewStartPosition(2);
        Assert.assertTrue(cache.cacheHit(parser.createGenomeLoc("1", 2, 100)), "Unexpected cache miss");
        Assert.assertEquals(cache.getCachedFeaturesUpToStopPosition(100), emptyRegion, "Should get back empty List for empty region");
    }

    /*********************************************************
     * End of direct testing on the FeatureCache inner class
     *********************************************************/
}
