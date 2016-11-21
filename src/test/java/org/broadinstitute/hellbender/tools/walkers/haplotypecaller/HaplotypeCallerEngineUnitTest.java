package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class HaplotypeCallerEngineUnitTest extends BaseTest {

    @Test
    public void testIsActive() throws IOException {
        final File testBam = new File(NA12878_20_21_WGS_bam);
        final File reference = new File(b37_reference_20_21);
        final SimpleInterval shardInterval = new SimpleInterval("20", 10000000, 10001000);
        final SimpleInterval paddedShardInterval = new SimpleInterval(shardInterval.getContig(), shardInterval.getStart() - 100, shardInterval.getEnd() + 100);
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

        // We expect isActive() to return 1.0 for the sites below, and 0.0 for all other sites
        final List<SimpleInterval> expectedActiveSites = Arrays.asList(
                new SimpleInterval("20", 9999996, 9999996),
                new SimpleInterval("20", 9999997, 9999997),
                new SimpleInterval("20", 10000117, 10000117),
                new SimpleInterval("20", 10000211, 10000211),
                new SimpleInterval("20", 10000439, 10000439),
                new SimpleInterval("20", 10000598, 10000598),
                new SimpleInterval("20", 10000694, 10000694),
                new SimpleInterval("20", 10000758, 10000758),
                new SimpleInterval("20", 10001019, 10001019)
        );

        try ( final ReadsDataSource reads = new ReadsDataSource(testBam.toPath());
              final ReferenceDataSource ref = new ReferenceFileSource(reference);
              final CachingIndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(reference);) {

            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, reads.getHeader(), referenceReader);

            List<ReadFilter> hcFilters = HaplotypeCallerEngine.makeStandardHCReadFilters();
            hcFilters.forEach(filter -> filter.setHeader(reads.getHeader()));
            ReadFilter hcCombinedFilter = hcFilters.get(0);
            for ( int i = 1; i < hcFilters.size(); ++i ) {
                hcCombinedFilter = hcCombinedFilter.and(hcFilters.get(i));
            }
            final Iterator<GATKRead> readIter = new ReadFilteringIterator(reads.query(paddedShardInterval), hcCombinedFilter);

            final LocusIteratorByState libs = new LocusIteratorByState(readIter, DownsamplingMethod.NONE, false, ReadUtils.getSamplesFromHeader(reads.getHeader()), reads.getHeader(), false);

            for ( final AlignmentContext pileup : libs ) {
                final SimpleInterval pileupInterval = new SimpleInterval(pileup.getLocation());
                final ReferenceContext pileupRefContext = new ReferenceContext(ref, pileupInterval);

                final ActivityProfileState isActiveResult = hcEngine.isActive(pileup, pileupRefContext, new FeatureContext(null, pileupInterval));

                final double expectedIsActiveValue = expectedActiveSites.contains(pileupInterval) ? 1.0 : 0.0;
                Assert.assertEquals(isActiveResult.isActiveProb(), expectedIsActiveValue, "Wrong isActive probability for site " + pileupInterval);
            }
        }
    }
}
