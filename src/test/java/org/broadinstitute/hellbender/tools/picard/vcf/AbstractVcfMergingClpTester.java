package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public abstract class AbstractVcfMergingClpTester extends CommandLineProgramTest {

	protected static final File TEST_DATA_PATH = new File(getTestDataDir(), "picard/vcf/trio");

    protected void runClp(final List<File> inputs, final File output) {
        runClp(inputs, output, Collections.<String>emptyList());
    }

    protected void runClp(final List<File> inputs, final File output, final List<String> otherArguments) {
        final List<String> argList = new ArrayList<>();
        for (final File input : inputs) {
            argList.add("--input");
            argList.add(input.getAbsolutePath());
        }
        argList.add("--output");
        argList.add(output.getAbsolutePath());
        argList.addAll(otherArguments);
        runCommandLine(argList);
    }

	@Test(expectedExceptions = {IllegalArgumentException.class, GATKException.class})
	public void testFailsOnDissimilarContigLists() {
		final File dissimilarContigs = new File(TEST_DATA_PATH, "CEUTrio-indels-dissimilar-contigs.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = new File("/dev/null/blah");
        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");

        runClp(Arrays.asList(dissimilarContigs, snpInputFile), output, indexing);
	}

	@Test(expectedExceptions = {IllegalArgumentException.class, GATKException.class})
	public void testFailsOnNoContigList() {
		final File contiglessIndelFile = new File(TEST_DATA_PATH, "CEUTrio-indels-no-contigs.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = new File("/dev/null");

        runClp(Arrays.asList(contiglessIndelFile, snpInputFile), output);
	}

    @Test(expectedExceptions = {IllegalArgumentException.class, GATKException.class})
	public void testFailsOnDissimilarSampleLists() {
		final File badSampleIndelFile = new File(TEST_DATA_PATH, "CEUTrio-indels-bad-samples.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = new File("/dev/null");

        runClp(Arrays.asList(badSampleIndelFile, snpInputFile), output);
	}

	@Test
	public void testMergeIndelsSnps() throws IOException {
		final File indelInputFile = new File(TEST_DATA_PATH, "CEUTrio-indels.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = BaseTest.createTempFile("merge-indels-snps-test-output.", ".vcf");
        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");

		final Queue<String> indelContigPositions = loadContigPositions(indelInputFile);
		final Queue<String> snpContigPositions = loadContigPositions(snpInputFile);

        runClp(Arrays.asList(indelInputFile, snpInputFile), output, indexing);
        validateSnpAndIndelResults(output, indelContigPositions, snpContigPositions);
    }

    /**
     * Make sure that the order of the output file is identical to the order
     * of the input files by iterating through the output, making sure that,
     * if the context is an indel (snp), the next genomic position in the indel
     * (snp) queue is the same. Also make sure that the context is in the order
     * specified by the input files.
     */
    private void validateSnpAndIndelResults(final File output, final Queue<String> indelContigPositions, final Queue<String> snpContigPositions) {
        final VCFFileReader outputReader = new VCFFileReader(output, false);
        final VariantContextComparator outputComparator = outputReader.getFileHeader().getVCFRecordComparator();
        VariantContext last = null;
        try (final CloseableIterator<VariantContext> iterator = outputReader.iterator()) {
            while (iterator.hasNext()) {
                final VariantContext outputContext = iterator.next();
                if (outputContext.isIndel())
                    Assert.assertEquals(getContigPosition(outputContext), indelContigPositions.poll());
                if (outputContext.isSNP())
                    Assert.assertEquals(getContigPosition(outputContext), snpContigPositions.poll());
                if (last != null) Assert.assertTrue(outputComparator.compare(last, outputContext) <= 0);
                last = outputContext;
            }
        }

        // We should have polled everything off the indel (snp) queues
        Assert.assertEquals(indelContigPositions.size(), 0);
        Assert.assertEquals(snpContigPositions.size(), 0);
    }

    @Test
	public void testMergeRandomScatter() throws IOException {
		final File zero = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-0.vcf");
		final File one = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-1.vcf");
		final File two = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-2.vcf");
		final File three = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-3.vcf");
		final File four = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-4.vcf");
		final File five = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-5.vcf");
        final List<File> inputs = Arrays.asList(zero, one, two, three, four, five);

		final List<Queue<String>> positionQueues = new ArrayList<>(6);
		positionQueues.add(0, loadContigPositions(zero));
		positionQueues.add(1, loadContigPositions(one));
		positionQueues.add(2, loadContigPositions(two));
		positionQueues.add(3, loadContigPositions(three));
		positionQueues.add(4, loadContigPositions(four));
		positionQueues.add(5, loadContigPositions(five));

        final List<String> indexing = Arrays.asList("--CREATE_INDEX", "false");
        final File output = BaseTest.createTempFile("random-scatter-test-output.", ".vcf");

        runClp(inputs, output, indexing);
        validateResultsForMultipleInputs(output, positionQueues);
	}

    private void validateResultsForMultipleInputs(final File output, final List<Queue<String>> positionQueues) {
        final VCFFileReader outputReader = new VCFFileReader(output, false);
        final VariantContextComparator outputComparator = outputReader.getFileHeader().getVCFRecordComparator();
        VariantContext last = null;
        try (final CloseableIterator<VariantContext> iterator = outputReader.iterator()) {
            while (iterator.hasNext()) {
                final VariantContext outputContext = iterator.next();
                final String position = getContigPosition(outputContext);
                for (final Queue<String> positionQueue : positionQueues) {
                    if (position.equals(positionQueue.peek())) {
                        positionQueue.poll();
                        break;
                    }
                }
                if (last != null) Assert.assertTrue(outputComparator.compare(last, outputContext) <= 0);
                last = outputContext;
            }
        }

        for (final Queue<String> positionQueue : positionQueues) {
            Assert.assertEquals(positionQueue.size(), 0);
        }
    }

    static Queue<String> loadContigPositions(final File inputFile) {
        final Queue<String> contigPositions = new LinkedList<>();
		try (final VCFFileReader reader = new VCFFileReader(inputFile, false);
             final CloseableIterator<VariantContext> iterator = reader.iterator()) {
                while (iterator.hasNext()) contigPositions.add(getContigPosition(iterator.next()));
        }
        return contigPositions;
	}

	static String getContigPosition(final VariantContext context) {
		return context.getContig() + "-" + Integer.toString(context.getStart());
	}
}
