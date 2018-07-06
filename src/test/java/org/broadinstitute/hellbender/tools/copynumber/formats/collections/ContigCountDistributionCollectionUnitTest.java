package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class ContigCountDistributionCollectionUnitTest {
    @Test
    public void test() {
        final String sampleName = "test-sample";
        final SimpleCountCollection counts = new SimpleCountCollection(
                new SimpleSampleLocatableMetadata(
                        sampleName,
                        new SAMSequenceDictionary(Arrays.asList(
                                new SAMSequenceRecord("chr1", 10000),
                                new SAMSequenceRecord("chr2", 10000)
                        ))),
                Arrays.asList(
                        new SimpleCount(new SimpleInterval("chr1", 1, 1000), 0),
                        new SimpleCount(new SimpleInterval("chr1", 1001, 2000), 1),
                        new SimpleCount(new SimpleInterval("chr1", 2001, 3000), 2),
                        new SimpleCount(new SimpleInterval("chr1", 3001, 4000), 3),
                        new SimpleCount(new SimpleInterval("chr1", 4001, 5000), 4),
                        new SimpleCount(new SimpleInterval("chr1", 5001, 6000), 0),
                        new SimpleCount(new SimpleInterval("chr1", 6001, 7000), 1),
                        new SimpleCount(new SimpleInterval("chr1", 7001, 8000), 2),
                        new SimpleCount(new SimpleInterval("chr1", 8001, 9000), 3),
                        new SimpleCount(new SimpleInterval("chr1", 9001, 10000), 10),
                        new SimpleCount(new SimpleInterval("chr2", 1, 1000), 0),
                        new SimpleCount(new SimpleInterval("chr2", 1001, 2000), 1),
                        new SimpleCount(new SimpleInterval("chr2", 2001, 3000), 2),
                        new SimpleCount(new SimpleInterval("chr2", 3001, 4000), 3),
                        new SimpleCount(new SimpleInterval("chr2", 4001, 5000), 4),
                        new SimpleCount(new SimpleInterval("chr2", 5001, 6000), 0),
                        new SimpleCount(new SimpleInterval("chr2", 6001, 7000), 1),
                        new SimpleCount(new SimpleInterval("chr2", 7001, 8000), 2),
                        new SimpleCount(new SimpleInterval("chr2", 8001, 9000), 3),
                        new SimpleCount(new SimpleInterval("chr2", 9001, 10000), 10)
                )
        );

        final Set<SimpleInterval> intervalSubset = new HashSet<>(counts.getIntervals());
        final int maximumCount = 5;
        final ContigCountDistributionCollection result = new ContigCountDistributionCollection(
                counts, intervalSubset, maximumCount);

        final Map<Integer, Integer> countDistributionExpected = new LinkedHashMap<>();
        countDistributionExpected.put(0, 2);
        countDistributionExpected.put(1, 2);
        countDistributionExpected.put(2, 2);
        countDistributionExpected.put(3, 2);
        countDistributionExpected.put(4, 1);
        countDistributionExpected.put(5, 0);

        Assert.assertEquals(result.getMetadata(), counts.getMetadata());
        Assert.assertEquals(result.getMaximumCount(), maximumCount);
        result.getRecords().forEach(
                ccd -> Assert.assertEquals(ccd.getCountDistribution(), countDistributionExpected));

        final File outputFile = IOUtils.createTempFile("contig-count-distribution", "tsv");
        result.write(outputFile);
        Assert.assertTrue(outputFile.exists());
    }
}