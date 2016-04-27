package org.broadinstitute.hellbender.tools.genome;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.datasources.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.exome.*;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SparkGenomeReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/genome");
    private static final File BAM_FILE = new File(TEST_FILE_DIR, "HCC1143_chr3_1K_11K.tiny.bam");
    private static final File REFERENCE_FILE = new File("src/test/resources/hg19mini.fasta");

    @Test
    public void testSparkGenomeReadCounts() {
        final File outputFile = createTempFile(BAM_FILE.getName(),".cov");
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BINSIZE_SHORT_NAME, "10000",
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        List<TargetCoverage> targetCoverages = loadAsTargetCoverages(outputFile);
        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 8);
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 16000);
        Assert.assertEquals(bedFeatureCollection.target(5).getName(), "target_3_10001_16000");
        Assert.assertEquals(targetCoverages.size(), bedFeatureCollection.targetCount());
    }

    @Test
    public void testSparkGenomeReadCountsBigBins() {
        final File outputFile = createTempFile(BAM_FILE.getName(), ".cov");
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BINSIZE_SHORT_NAME, "16000",
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        List<TargetCoverage> targetCoverages = loadAsTargetCoverages(outputFile);
        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 4);
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 16000);
        Assert.assertEquals(bedFeatureCollection.target(2).getName(), "target_3_1_16000");
        Assert.assertEquals(targetCoverages.size(), bedFeatureCollection.targetCount());

    }

    @Test
    public void testSparkGenomeReadCountsSmallBins() {
        final File outputFile = createTempFile(BAM_FILE.getName(), ".cov");
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BINSIZE_SHORT_NAME, "2000",
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);

        // Proportional Coverage
        final List<TargetCoverage> targetProportionalCoverages = loadAsTargetCoverages(outputFile);
        Assert.assertTrue(targetProportionalCoverages.stream().anyMatch(t -> Math.abs(t.getCoverage()) > 1e-10));

        // The reads are all in three bins of contig 3 with values {.5, .25, .25}
        Assert.assertTrue(targetProportionalCoverages.stream().filter(t -> t.getContig().equals("3")).anyMatch(t -> Math.abs(t.getCoverage()) > .2));
        Assert.assertTrue(Math.abs(targetProportionalCoverages.stream().filter(t -> t.getContig().equals("3")).mapToDouble(t -> t.getCoverage()).sum() - 1.0) < 1e-10);

        // raw coverage
        final List<TargetCoverage> targetCoverages = loadAsTargetCoverages(new File(outputFile.getAbsolutePath() + SparkGenomeReadCounts.RAW_COV_OUTPUT_EXTENSION));

        Assert.assertTrue(targetCoverages.stream().anyMatch(t -> Math.abs(t.getCoverage()) > 1e-10));

        // The reads are all in three bins of contig 3 with values
        Assert.assertEquals(targetCoverages.stream().filter(t -> t.getContig().equals("3")).filter(t -> Math.abs(t.getCoverage()) >= 1).count(), 3);

        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        final TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 16000/2000 * 4); // 4 is the number of contigs in the fasta file
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 4000);
        Assert.assertEquals(bedFeatureCollection.target(2).getName(), "target_1_4001_6000");
        Assert.assertEquals(bedFeatureCollection.target(8).getName(), "target_2_1_2000");
        Assert.assertEquals(bedFeatureCollection.target(17).getName(), "target_3_2001_4000");
        Assert.assertEquals(targetProportionalCoverages.size(), bedFeatureCollection.targetCount());
    }

    private List<TargetCoverage> loadAsTargetCoverages(File outputFile) {
        ReadCountCollection targetProportionalCoveragesAsRCC = null;
        try {
            targetProportionalCoveragesAsRCC = ReadCountCollectionUtils.parse(outputFile);
        } catch (final IOException ioe) {
            Assert.fail("IO Exception in automated test.  Possible misconfiguration?", ioe);
        }

        // Convert to target coverages for convenience in testing.
        return ReadCountCollectionUtils.convertToTargetCoverageList(targetProportionalCoveragesAsRCC);
    }

    @Test
    public void testDropContigsFromSequence() throws IOException{
        final Set<String> testContigsToDrop = new HashSet<>(Arrays.asList("3"));
        ReferenceFileSource testFasta = new ReferenceFileSource(REFERENCE_FILE.getAbsolutePath());
        final SAMSequenceDictionary originalSequenceDictionary = testFasta.getReferenceSequenceDictionary(null);
        SAMSequenceDictionary finalSequenceDictionary = SparkGenomeReadCounts.dropContigsFromSequence(originalSequenceDictionary, testContigsToDrop);
        Assert.assertEquals(finalSequenceDictionary.size(), 3);

    }
}
