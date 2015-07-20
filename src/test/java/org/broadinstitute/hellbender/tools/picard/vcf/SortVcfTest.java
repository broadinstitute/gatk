package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Test class for SortVCF. Several tests lso conducted by AbstractVcfMergingClpTester
 *
 * Created by bradt on 9/3/14.
 */
public final class SortVcfTest extends AbstractVcfMergingClpTester {

    @Test
    public void testPresortedFile() throws IOException {
        final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = BaseTest.createTempFile("sort-presorted-test-output.", ".vcf");
        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");

        final int numberOfVariantContexts = loadContigPositions(snpInputFile).size();

        runClp(Arrays.asList(snpInputFile), output, indexing);
        validateSortingResults(output, numberOfVariantContexts);
    }

    @Test
    public void testSingleScrambledFile() throws IOException {
        final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps-scrambled.1.vcf");
        final File output = BaseTest.createTempFile("sort-single-scrambled-test-output.", ".vcf");
        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");

        final int numberOfVariantContexts = loadContigPositions(snpInputFile).size();

        runClp(Arrays.asList(snpInputFile), output, indexing);
        validateSortingResults(output, numberOfVariantContexts);
    }

    @Test
    public void testTwoScrambledSnpFiles() throws IOException {
        final File inputFile1 = new File(TEST_DATA_PATH, "CEUTrio-snps-scrambled.1.vcf");
        final File inputFile2 = new File(TEST_DATA_PATH, "CEUTrio-snps-scrambled.2.vcf");
        final File output = BaseTest.createTempFile("sort-multiple-scrambled-test-output.", ".vcf");
        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");

        final int numberOfVariantContexts = loadContigPositions(inputFile1).size() + loadContigPositions(inputFile2).size();

        runClp(Arrays.asList(inputFile1, inputFile2), output, indexing);
        validateSortingResults(output, numberOfVariantContexts);
    }

    @Test
    public void testScrambledSnpsAndOrderedIndels() throws IOException {
        final File indelInputFile = new File(TEST_DATA_PATH, "CEUTrio-indels.vcf");
        final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps-scrambled.1.vcf");
        final File output = BaseTest.createTempFile("sort-scrambled-indels-snps-test-output.", ".vcf");
        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");

        final int numberOfVariantContexts = loadContigPositions(indelInputFile).size() + loadContigPositions(snpInputFile).size();

        runClp(Arrays.asList(indelInputFile, snpInputFile), output, indexing);
        validateSortingResults(output, numberOfVariantContexts);
    }

    @Test
    public void testScrambledSnpsAndScrambledIndels() throws IOException {
        final File indelInputFile = new File(TEST_DATA_PATH, "CEUTrio-indels-scrambled.1.vcf");
        final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps-scrambled.1.vcf");
        final File output = BaseTest.createTempFile("merge-indels-snps-test-output.", ".vcf");
        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");;

        final int numberOfVariantContexts = loadContigPositions(indelInputFile).size() + loadContigPositions(snpInputFile).size();

        runClp(Arrays.asList(indelInputFile, snpInputFile), output, indexing);
        validateSortingResults(output, numberOfVariantContexts);
    }


    /**
     * Checks the ordering and total number of variant context entries in the specified output VCF file.
     * Does NOT check explicitly that the VC genomic positions match exactly those from the inputs. We assume this behavior from other tests.
     *
     * @param output VCF file representing the output of SortVCF
     * @param expectedVariantContextCount the total number of variant context entries from all input files that were merged/sorted
     */
    private void validateSortingResults(final File output, final int expectedVariantContextCount) {
        final VCFFileReader outputReader = new VCFFileReader(output, false);
        final VariantContextComparator outputComparator = outputReader.getFileHeader().getVCFRecordComparator();
        VariantContext last = null;
        int variantContextCount = 0;
        try (final CloseableIterator<VariantContext> iterator = outputReader.iterator()) {
            while (iterator.hasNext()) {
                final VariantContext outputContext = iterator.next();
                if (last != null) Assert.assertTrue(outputComparator.compare(last, outputContext) <= 0);
                last = outputContext;
                variantContextCount++;
            }
        }
        Assert.assertEquals(variantContextCount, expectedVariantContextCount);
    }
}
