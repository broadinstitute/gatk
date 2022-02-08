package org.broadinstitute.hellbender.engine.cache;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class LocatableCacheUnitTest {

    /********************************************************************************************
     * LocatableCache tests using DrivingFeatureInputCacheStrategy to test basic cache operations
     ********************************************************************************************/

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
        //suppressing deprecation; function required because it's part of the implemented class
        @Override
        @SuppressWarnings("deprecation")
        @Deprecated
        public String getChr() {
            return getContig();
        }

        @Override
        public String getContig() {
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
            return chr.equals(otherFeature.getContig()) && start == otherFeature.getStart() && end == otherFeature.getEnd();
        }

        @Override
        public String toString() {
            return chr + ":" + start + "-" + end;   // (to improve output on test failures involving this class)
        }
    }

    private LocatableCache<ArtificialTestFeature> initializeFeatureCache(final List<ArtificialTestFeature> features, final String cacheContig, final int cacheStart, final int cacheEnd, final int cacheLookAhead ) {
        LocatableCache<ArtificialTestFeature> cache = new LocatableCache<>(
                "test",
                new DrivingFeatureInputCacheStrategy<>(
                        cacheLookAhead,
                        (Locatable newCacheInterval) -> features.iterator())
        );

        cache.queryAndPrefetch(new SimpleInterval(cacheContig, cacheStart, cacheEnd));
        return cache;
    }

    @DataProvider(name = "FeatureCacheInitialQueryDataProvider")
    public Object[][] getFeatureCacheInitialQueryData() {
        return new Object[][] {
                { Arrays.asList(
                        new ArtificialTestFeature("1", 1, 100),
                        new ArtificialTestFeature("1", 50, 150),
                        new ArtificialTestFeature("1", 200, 300),
                        new ArtificialTestFeature("1", 350, 400)),
                        "1", 1, 400 },
                { Arrays.asList(new ArtificialTestFeature("1", 1, 100)), "1", 1, 100 },
                { Collections.<ArtificialTestFeature>emptyList(), "1", 1, 1 }
        };
    }

    @Test(dataProvider = "FeatureCacheInitialQueryDataProvider")
    public void testFeatureCacheInitialQuery( final List<ArtificialTestFeature> features, final String cacheContig, final int cacheStart, final int cacheEnd) {
        LocatableCache<ArtificialTestFeature> cache = initializeFeatureCache(features, cacheContig, cacheStart, cacheEnd, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);

        List<ArtificialTestFeature> cachedFeatures = cache.getCachedLocatables(new SimpleInterval(cache.getContig(), cache.getCacheStart(), cacheEnd));
        Assert.assertEquals(cache.getContig(), cacheContig, "Wrong contig reported by cache after fill");
        Assert.assertEquals(cache.getCacheStart(), cacheStart, "Wrong start position reported by cache after fill");
        Assert.assertEquals(cache.getCacheEnd(), cacheEnd + FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, "Wrong stop position reported by cache after fill");
        Assert.assertEquals(cachedFeatures, features, "Wrong Features in cache after fill()");
    }

    @DataProvider(name = "FeatureCacheHitDetectionDataProvider")
    public Object[][] getFeatureCacheHitDetectionData() {
        List<ArtificialTestFeature> features = Arrays.asList(new ArtificialTestFeature("1", 1, 100),
                new ArtificialTestFeature("1", 50, 150),
                new ArtificialTestFeature("1", 200, 300));
        // initialize cache with 0 lookahead bases
        LocatableCache<ArtificialTestFeature> cache = initializeFeatureCache(features, "1", 50, 250, 0);

        return new Object[][] {
                // Exact match for cache boundaries
                { cache, new SimpleInterval("1", 50, 250), true },
                // Interval completely contained within cache boundaries
                { cache, new SimpleInterval("1", 100, 200), true },
                // Interval left-aligned with cache boundaries
                { cache, new SimpleInterval("1", 50, 100), true },
                // Interval right-aligned with cache boundaries
                { cache, new SimpleInterval("1", 200, 250), true },
                // Interval overlaps, but is off the left edge of cache boundaries
                { cache, new SimpleInterval("1", 49, 100), false },
                // Interval overlaps, but is off the right edge of cache boundaries
                { cache, new SimpleInterval("1", 200, 251), false },
                // Interval does not overlap and is to the left of cache boundaries
                { cache, new SimpleInterval("1", 1, 40), false },
                // Interval does not overlap and is to the right of cache boundaries
                { cache, new SimpleInterval("1", 300, 350), false },
                // Interval is on different contig
                { cache, new SimpleInterval("2", 50, 250), false }
        };
    }

    @Test(dataProvider = "FeatureCacheHitDetectionDataProvider")
    public void testFeatureCacheHitDetection( final LocatableCache<ArtificialTestFeature> cache,
                                       final SimpleInterval testInterval, final boolean cacheHitExpectedResult ) {
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
        // initialize the feature cache with 0 lookahead bases
        LocatableCache<ArtificialTestFeature> cache = initializeFeatureCache(feats, "1", 1, 200, 0);

        // Pairing of start position to which to trimCache the cache with the List of Features we expect to see
        // in the cache after trimming
        List<Pair<Integer, List<ArtificialTestFeature>>> trimOperations = Arrays.asList(
                Pair.of(1, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(2, Arrays.asList(feats.get(1), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(3, Arrays.asList(feats.get(1), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(4, Arrays.asList(feats.get(1), feats.get(3), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(5, Arrays.asList(feats.get(1), feats.get(3), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(6, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(10, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(8), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(11, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(50, Arrays.asList(feats.get(1), feats.get(3), feats.get(7), feats.get(9), feats.get(10), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(51, Arrays.asList(feats.get(1), feats.get(9), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(100, Arrays.asList(feats.get(1), feats.get(9), feats.get(11), feats.get(12), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(101, Arrays.asList(feats.get(11), feats.get(13), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(111, Arrays.asList(feats.get(11), feats.get(14), feats.get(15), feats.get(16))),
                Pair.of(151, Arrays.asList(feats.get(11), feats.get(14), feats.get(16))),
                Pair.of(151, Arrays.asList(feats.get(11), feats.get(14), feats.get(16))),
                Pair.of(200, Arrays.asList(feats.get(11), feats.get(14)))
        );

        return new Object[][] {
                { cache, trimOperations }
        };
    }

    @Test(dataProvider = "FeatureCacheTrimmingDataProvider")
    public void testFeatureCacheTrimming( final LocatableCache<ArtificialTestFeature> cache, final List<Pair<Integer, List<ArtificialTestFeature>>> trimOperations ) {
        // Repeatedly trimCache the cache to ever-increasing start positions, and verify after each trimCache operation
        // that the cache holds the correct Features in the correct order
        for ( Pair<Integer, List<ArtificialTestFeature>> trimOperation : trimOperations ) {
            final int trimPosition = trimOperation.getLeft();
            final List<ArtificialTestFeature> expectedFeatures = trimOperation.getRight();

            cache.queryAndPrefetch(new SimpleInterval(cache.getContig(), trimPosition, cache.getCacheEnd()));

            final List<ArtificialTestFeature> actualFeatures = cache.getCachedLocatables(new SimpleInterval(cache.getContig(), cache.getCacheStart(), cache.getCacheEnd()));
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
        LocatableCache<ArtificialTestFeature> cache = initializeFeatureCache(feats, "1", 1, 100, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES);

        // Pairing of end position with which to bound cache retrieval with the List of Features we expect to see
        // after retrieval
        List<Pair<Integer, List<ArtificialTestFeature>>> retrievalOperations = Arrays.asList(
                Pair.of(100, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(80, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(79, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9))),
                Pair.of(80, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9), feats.get(10))),
                Pair.of(75, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8), feats.get(9))),
                Pair.of(74, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7), feats.get(8))),
                Pair.of(54, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7))),
                Pair.of(52, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6), feats.get(7))),
                Pair.of(51, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5), feats.get(6))),
                Pair.of(50, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4), feats.get(5))),
                Pair.of(49, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4))),
                Pair.of(10, Arrays.asList(feats.get(0), feats.get(1), feats.get(2), feats.get(3), feats.get(4))),
                Pair.of(9, Arrays.asList(feats.get(0), feats.get(1), feats.get(2))),
                Pair.of(5, Arrays.asList(feats.get(0), feats.get(1), feats.get(2))),
                Pair.of(4, Arrays.asList(feats.get(0), feats.get(1))),
                Pair.of(1, Arrays.asList(feats.get(0), feats.get(1)))
        );

        return new Object[][] {
                { cache, retrievalOperations }
        };
    }

    @Test(dataProvider = "FeatureCacheRetrievalDataProvider")
    public void testFeatureCacheFeatureRetrieval( final LocatableCache<ArtificialTestFeature> cache, final List<Pair<Integer, List<ArtificialTestFeature>>> retrievalOperations ) {
        for ( Pair<Integer, List<ArtificialTestFeature>> retrievalOperation: retrievalOperations ) {
            final int stopPosition = retrievalOperation.getLeft();
            final List<ArtificialTestFeature> expectedFeatures = retrievalOperation.getRight();

            final List<ArtificialTestFeature> actualFeatures = cache.getCachedLocatables(new SimpleInterval(cache.getContig(), cache.getCacheStart(), stopPosition));
            Assert.assertEquals(actualFeatures, expectedFeatures, "Wrong Features returned in retrieval operation with stop position " + stopPosition);
        }
    }

    /**
     * Test caching a region with no Features. This should work (we should avoid going to disk
     * to look for new records when querying within such a region).
     */
    @Test
    public void testHandleCachingOfEmptyRegion() {
        List<ArtificialTestFeature> emptyRegion = new ArrayList<>();
        // initialize cache with 0 lookahead bases
        LocatableCache<ArtificialTestFeature> cache = new LocatableCache<>(
                "test",
                new DrivingFeatureInputCacheStrategy<ArtificialTestFeature>(
                        0, (Locatable newInterval) -> emptyRegion.iterator()));

        // prime the cache with an empty iterator
        cache.queryAndPrefetch(new SimpleInterval("1", 1, 100));

        Assert.assertTrue(cache.isEmpty(), "Cache should be empty");
        Assert.assertTrue(cache.cacheHit(new SimpleInterval("1", 1, 100)), "Unexpected cache miss");
        Assert.assertTrue(cache.cacheHit(new SimpleInterval("1", 2, 99)), "Unexpected cache miss");

        // query the empty region to force a (no-op) trimCache, after which we still expect a cacheHit
        Assert.assertEquals(
                cache.queryAndPrefetch(
                    new SimpleInterval(cache.getContig(), 2, cache.getCacheEnd())),
                    emptyRegion);
        Assert.assertTrue(cache.cacheHit(new SimpleInterval("1", 2, 100)), "Unexpected cache miss");

        // now query again and make sure we don't do another query
        Assert.assertEquals(
                cache.queryAndPrefetch(
                        new SimpleInterval(cache.getContig(), cache.getCacheStart(), 100)),
                emptyRegion,
                "Should get back empty List for empty region");
    }

    /***************************************************************************************************
     * End of LocatableCache tests using DrivingFeatureInputCacheStrategy to test basic cache operations
     ***************************************************************************************************/

    /********************************************************
     * LocatableCache tests using SideReadInputCacheStrategy
     ********************************************************/

    @DataProvider(name = "ReadsWithUnmappedMates")
    public Object[][] getReadsWithUnmappedMates() {
        // hit every code path in SideReadInputCacheStrategy:
        //
        // mate pair appears at an early locus/late locus (i.e., before/after other reads)
        // mate pairs with mapped first/unmapped first
        // query that trims/doesn't trim the pair from the cache

        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();

        final SimpleInterval earlyInterval = new SimpleInterval("1", 1, 10);
        final SimpleInterval lateInterval = new SimpleInterval("1", 20, 30);

        final GATKRead earlyMappedMate = ArtificialReadUtils.createArtificialRead(samHeader, "earlyMatePair", 0, earlyInterval.getStart(), 10);
        earlyMappedMate.setIsPaired(true);
        final GATKRead earlyUnmappedMate = ArtificialReadUtils.createArtificialRead(samHeader, "earlyMatePair", 0, earlyInterval.getStart(), 10);
        earlyUnmappedMate.setIsPaired(true);
        earlyUnmappedMate.setIsUnmapped();

        final GATKRead lateMappedMate = ArtificialReadUtils.createArtificialRead(samHeader, "lateMatePair", 0, lateInterval.getStart(), 10);
        lateMappedMate.setIsPaired(true);
        final GATKRead lateUnmappedMate = ArtificialReadUtils.createArtificialRead(samHeader, "lateMatePair", 0, lateInterval.getStart(), 10);
        lateUnmappedMate.setIsPaired(true);
        lateUnmappedMate.setIsUnmapped();

        final GATKRead earlyOtherRead1 = ArtificialReadUtils.createArtificialRead(samHeader, "earlyOtherRead1", 0, earlyInterval.getStart(), 10);
        final GATKRead earlyOtherRead2 = ArtificialReadUtils.createArtificialRead(samHeader, "earlyOtherRead2", 0, earlyInterval.getStart(), 10);
        final GATKRead lateOtherRead1 = ArtificialReadUtils.createArtificialRead(samHeader, "lateOtherRead1", 0, lateInterval.getStart(), 10);
        final GATKRead lateOtherRead2 = ArtificialReadUtils.createArtificialRead(samHeader, "lateOtherRead2", 0, lateInterval.getStart(), 10);

        // add in coord sorted order:
        // unmapped mate, followed by mapped mate, followed by other reads
        final List<GATKRead> earlyUnmappedMateFirst = new ArrayList<>();
        earlyUnmappedMateFirst.add(earlyUnmappedMate);      // unmapped mate
        earlyUnmappedMateFirst.add(earlyMappedMate);
        earlyUnmappedMateFirst.add(lateOtherRead1);
        earlyUnmappedMateFirst.add(lateOtherRead2);

        // mapped mate, followed by unmapped mate, followed by other reads
        final List<GATKRead> earlyUnmappedMateSecond = new ArrayList<>();
        earlyUnmappedMateSecond.add(earlyMappedMate);
        earlyUnmappedMateSecond.add(earlyUnmappedMate);     // unmapped mate
        earlyUnmappedMateSecond.add(lateOtherRead1);
        earlyUnmappedMateSecond.add(lateOtherRead2);

        // mate pair appears after other reads, unmapped mate first
        final List<GATKRead> lateUnmappedMateFirst = new ArrayList<>();
        lateUnmappedMateFirst.add(earlyOtherRead1);
        lateUnmappedMateFirst.add(earlyOtherRead2);
        lateUnmappedMateFirst.add(lateUnmappedMate);         // unmapped mate
        lateUnmappedMateFirst.add(lateMappedMate);

        // mate pair appears after other reads, mapped mate first
        final List<GATKRead> lateUnmappedMateSecond = new ArrayList<>();
        lateUnmappedMateSecond.add(earlyOtherRead1);
        lateUnmappedMateSecond.add(earlyOtherRead2);
        lateUnmappedMateSecond.add(lateMappedMate);
        lateUnmappedMateSecond.add(lateUnmappedMate);         // unmapped mate

        return new Object[][]{
                // cache contents, lookaheadBases, initial cache triggering interval,
                //      query interval, triggering interval result size, index of unmapped, query result size, empty at end

                // issue a query past the end of the mate pair to force them both to be trimmed from the cache
                { earlyUnmappedMateFirst, 50, earlyInterval, new SimpleInterval("1", 20, 50), 2, 0, 2, false},

                // issue a query that overlaps both reads in the pair to cause them both to remain in the cache
                { earlyUnmappedMateSecond, 50, earlyInterval, new SimpleInterval("1", 1, 50), 2, 1, 4, false },

                // issue a query that overlaps only the MAPPED mate in the pair to cause them both  to remain in the cache
                { earlyUnmappedMateSecond, 50, earlyInterval, new SimpleInterval("1", 8, 50), 2, 1, 3, false },

                // issue a query that includes the mate pair to force them both to be trimmed from the cache
                { earlyUnmappedMateSecond, 50, earlyInterval, new SimpleInterval("1", 30, 40), 2, 1, 0, true },

                // issue a query past the end of mate pair to force them both to be trimmed from the cache
                { lateUnmappedMateFirst, 50, earlyInterval, new SimpleInterval("1", 40, 50), 2, -1, 0, true },

                // issue a query that overlaps only the MAPPED mate, so both remain in the cache
                { lateUnmappedMateSecond, 50, earlyInterval, new SimpleInterval("1", 25, 50), 2, -1, 1, false }
        };
    }

    @Test(dataProvider = "ReadsWithUnmappedMates")
    public void testReadCacheTrimmingWithUnmappedMates(
            final List<GATKRead> reads,
            final int lookAheadBases,
            final SimpleInterval triggeringQueryInterval,
            final SimpleInterval queryInterval,
            final int returnedByTriggeringInterval,
            final int indexOfUnmapped,
            final int returnedByQueryIntervalQuery,
            final boolean expectedEmpty
    ) {
        final LocatableCache<GATKRead> cache =
                new LocatableCache<>(
                    "test",
                    new SideReadInputCacheStrategy<>(
                            lookAheadBases,
                            (Locatable interval) -> reads.iterator()
                    )
                );

        List<GATKRead> resultReads = cache.queryAndPrefetch(triggeringQueryInterval);
        Assert.assertEquals(resultReads.size(), returnedByTriggeringInterval);
        Assert.assertTrue(
                indexOfUnmapped == -1 ||
                    (resultReads.get(0).isPaired() && resultReads.get(0).isUnmapped())
                        || (resultReads.get(0).isPaired() && resultReads.get(1).mateIsUnmapped()));

        // issue the query and see how many records are returned
        resultReads = cache.queryAndPrefetch(queryInterval);
        Assert.assertEquals(resultReads.size(), returnedByQueryIntervalQuery);

        Assert.assertEquals(cache.isEmpty(), expectedEmpty);
    }

}
