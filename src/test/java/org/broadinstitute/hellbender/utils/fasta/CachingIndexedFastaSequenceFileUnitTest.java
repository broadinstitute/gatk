package org.broadinstitute.hellbender.utils.fasta;


// the imports for unit testing.


import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.Level;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.BaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Basic unit test for CachingIndexedFastaSequenceFile
 */
public final class CachingIndexedFastaSequenceFileUnitTest extends GATKBaseTest {
    private Path simpleFasta = Paths.get(publicTestDir + "/exampleFASTA.fasta");
    private Path simpleBGZippedFasta = Paths.get(publicTestDir + "/exampleFASTA.fasta.gz");
    private static final int STEP_SIZE = 1;

    private static final List<Integer> QUERY_SIZES = Arrays.asList(1, 10, 100);
    private static final List<Integer> CACHE_SIZES = Arrays.asList(-1, 100, 1000);

    @DataProvider(name = "fastas")
    public Object[][] createData1() {
        List<Object[]> params = new ArrayList<>();
        // the tests compare to {@link simpleFasta} as the uncached version of the file
        // you can't add arbitrary fastas here without refactoring
        for ( Path fasta : Arrays.asList(simpleFasta, simpleBGZippedFasta) ) {
            for ( int cacheSize : CACHE_SIZES ) {
                for ( int querySize : QUERY_SIZES ) {
                    params.add(new Object[]{fasta, simpleFasta, cacheSize, querySize});
                }
            }
        }

        return params.toArray(new Object[][]{});
    }

    private static long getCacheSize(final long cacheSizeRequested) {
        return cacheSizeRequested == -1 ? CachingIndexedFastaSequenceFile.DEFAULT_CACHE_SIZE : cacheSizeRequested;
    }

    @Test(dataProvider = "fastas")
    public void testCachingIndexedFastaReaderSequential1(Path fasta, Path unzipped, int cacheSize, int querySize) throws IOException {
        try(final CachingIndexedFastaSequenceFile caching = new CachingIndexedFastaSequenceFile(fasta, getCacheSize(cacheSize), true, false)) {

            SAMSequenceRecord contig = caching.getSequenceDictionary().getSequence(0);
            logger.debug(String.format("Checking contig %s length %d with cache size %d and query size %d",
                                      contig.getSequenceName(), contig.getSequenceLength(), cacheSize, querySize));
            testSequential(caching, unzipped, querySize);
        }
    }

    private void testSequential(final CachingIndexedFastaSequenceFile caching, final Path unzipped, final int querySize) throws IOException {
        Assert.assertTrue(caching.isPreservingCase(), "testSequential only works for case preserving CachingIndexedFastaSequenceFile readers");

        //use an unzipped version of the fasta because it's slow to repeatedly query the zipped version without caching
        try( final ReferenceSequenceFile uncached = ReferenceSequenceFileFactory.getReferenceSequenceFile(unzipped)) {

            SAMSequenceRecord contig = uncached.getSequenceDictionary().getSequence(0);
            for (int i = 0; i < contig.getSequenceLength(); i += STEP_SIZE) {
                int start = i;
                int stop = start + querySize;
                if (stop <= contig.getSequenceLength()) {
                    ReferenceSequence cachedVal = caching.getSubsequenceAt(contig.getSequenceName(), start, stop);
                    ReferenceSequence uncachedVal = uncached.getSubsequenceAt(contig.getSequenceName(), start, stop);

                    Assert.assertEquals(cachedVal.getName(), uncachedVal.getName());
                    Assert.assertEquals(cachedVal.getContigIndex(), uncachedVal.getContigIndex());
                    Assert.assertEquals(cachedVal.getBases(), uncachedVal.getBases());
                }
            }

            // asserts for efficiency.  We are going to make contig.length / STEP_SIZE queries
            // at each of range: start -> start + querySize against a cache with size of X.
            // we expect to hit the cache each time range falls within X.  We expect a hit
            // on the cache if range is within X.  Which should happen at least (X - query_size * 2) / STEP_SIZE
            // times.
            final int minExpectedHits = (int) Math.floor(
                    (Math.min(caching.getCacheSize(), contig.getSequenceLength()) - querySize * 2.0) / STEP_SIZE);
            caching.printEfficiency(Level.DEBUG);
            Assert.assertTrue(caching.getCacheHits() >= minExpectedHits,
                              "Expected at least " + minExpectedHits + " cache hits but only got " + caching.getCacheHits());
        }

    }

    // Tests grabbing sequences around a middle cached value.
    @Test(dataProvider = "fastas")
    public void testCachingIndexedFastaReaderTwoStage(Path fasta, Path unzipped, int cacheSize, int querySize) throws IOException {
        try(final ReferenceSequenceFile uncached = ReferenceSequenceFileFactory.getReferenceSequenceFile(unzipped);
            final CachingIndexedFastaSequenceFile caching = new CachingIndexedFastaSequenceFile(fasta, getCacheSize(cacheSize), true, false)) {

            SAMSequenceRecord contig = uncached.getSequenceDictionary().getSequence(0);

            int middleStart = (contig.getSequenceLength() - querySize) / 2;
            int middleStop = middleStart + querySize;

            logger.debug(String.format(
                    "Checking contig %s length %d with cache size %d and query size %d with intermediate query",
                    contig.getSequenceName(), contig.getSequenceLength(), cacheSize, querySize));

            for (int i = 0; i < contig.getSequenceLength(); i += 10) {
                int start = i;
                int stop = start + querySize;
                if (stop <= contig.getSequenceLength()) {
                    ReferenceSequence grabMiddle = caching.getSubsequenceAt(contig.getSequenceName(), middleStart,
                                                                            middleStop);
                    ReferenceSequence cachedVal = caching.getSubsequenceAt(contig.getSequenceName(), start, stop);
                    ReferenceSequence uncachedVal = uncached.getSubsequenceAt(contig.getSequenceName(), start, stop);

                    Assert.assertEquals(cachedVal.getName(), uncachedVal.getName());
                    Assert.assertEquals(cachedVal.getContigIndex(), uncachedVal.getContigIndex());
                    Assert.assertEquals(cachedVal.getBases(), uncachedVal.getBases());
                }
            }
        }
    }

    // make sure some bases are lower case and some are upper case
    @Test
    public void testMixedCasesInExample() throws IOException {
        try(final IndexedFastaSequenceFile original = new IndexedFastaSequenceFile(new File(exampleFASTA));
            final CachingIndexedFastaSequenceFile casePreserving = new CachingIndexedFastaSequenceFile(IOUtils.getPath(exampleFASTA), true);
            final CachingIndexedFastaSequenceFile allUpper = new CachingIndexedFastaSequenceFile(IOUtils.getPath(exampleFASTA), CachingIndexedFastaSequenceFile.DEFAULT_CACHE_SIZE, false, true);
        ) {

            int nMixedCase = 0;
            for (SAMSequenceRecord contig : original.getSequenceDictionary().getSequences()) {
                nMixedCase += mixedCasesTestHelper(original, casePreserving, allUpper, contig.getSequenceName(), -1, -1);

                final int step = 100;
                for (int lastPos = step; lastPos < contig.getSequenceLength(); lastPos += step) {
                    mixedCasesTestHelper(original, casePreserving, allUpper, contig.getSequenceName(), lastPos - step, lastPos);
                }
            }


            Assert.assertTrue(nMixedCase > 0, "No mixed cases sequences found in file.  Unexpected test state");
        }
    }

    private int mixedCasesTestHelper(final ReferenceSequenceFile original,
                                     final ReferenceSequenceFile casePreserving,
                                     final ReferenceSequenceFile allUpper,
                                     final String contig, final int start, final int stop ) {
        final String orig = fetchBaseString(original, contig, start, stop);
        final String keptCase = fetchBaseString(casePreserving, contig, start, stop);
        final String upperCase = fetchBaseString(allUpper, contig, start, stop);

        final String origToUpper = orig.toUpperCase();
        if ( ! orig.equals(origToUpper) ) {
            Assert.assertEquals(keptCase, orig, "Case preserving operation not equal to the original case for contig " + contig);
            Assert.assertEquals(upperCase, origToUpper, "All upper case reader not equal to the uppercase of original case for contig " + contig);
            return 1;
        } else {
            return 0;
        }
    }

    private String fetchBaseString(final ReferenceSequenceFile reader, final String contig, final int start, final int stop) {
        if ( start == -1 ) {
            return new String(reader.getSequence(contig).getBases());
        } else {
            return new String(reader.getSubsequenceAt(contig, start, stop).getBases());
        }
    }

    @Test
    public void testIupacChanges() throws IOException {
        final String testFasta = publicTestDir + "iupacFASTA.fasta";
        try(final CachingIndexedFastaSequenceFile iupacPreserving = new CachingIndexedFastaSequenceFile(IOUtils.getPath(testFasta), CachingIndexedFastaSequenceFile.DEFAULT_CACHE_SIZE, false, true);
            final CachingIndexedFastaSequenceFile makeNs = new CachingIndexedFastaSequenceFile(IOUtils.getPath(testFasta), CachingIndexedFastaSequenceFile.DEFAULT_CACHE_SIZE, true, false)) {

            int preservingNs = 0;
            int changingNs = 0;
            for (SAMSequenceRecord contig : iupacPreserving.getSequenceDictionary().getSequences()) {
                final String sPreserving = fetchBaseString(iupacPreserving, contig.getSequenceName(), 0, 15000);
                preservingNs += StringUtils.countMatches(sPreserving, "N");

                final String sChanging = fetchBaseString(makeNs, contig.getSequenceName(), 0, 15000);
                changingNs += StringUtils.countMatches(sChanging, "N");
            }

            Assert.assertEquals(changingNs, preservingNs + 4);
        }
    }

    @Test(expectedExceptions = {UserException.class})
    public void testFailOnBadBase() throws IOException {
        final String testFasta = publicTestDir + "problematicFASTA.fasta";
        try(final CachingIndexedFastaSequenceFile fasta = new CachingIndexedFastaSequenceFile(IOUtils.getPath(testFasta))) {

            for (SAMSequenceRecord contig : fasta.getSequenceDictionary().getSequences()) {
                fetchBaseString(fasta, contig.getSequenceName(), -1, -1);
            }
        }
    }

    @Test(expectedExceptions = {UserException.MissingReference.class})
    public void testNonexistentReference() {
        new CachingIndexedFastaSequenceFile(BaseTest.getSafeNonExistentPath("NonexistentReference.fasta"));
    }

}
