package org.broadinstitute.hellbender.tools.genome;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.TargetCollectionUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SparkGenomeReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/genome");
    private static final File BAM_FILE = new File(TEST_FILE_DIR, "HCC1143_chr3_1K_11K.tiny.bam");
    private static final File REFERENCE_FILE = new File("src/test/resources/hg19mini.fasta");

    @Test
    public void testSparkGenomeReadCounts() throws IOException {
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
        final ReadCountCollection coverage = ReadCountCollectionUtils.parse(outputFile);
        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 8);
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 16000);
        Assert.assertEquals(bedFeatureCollection.target(5).getName(), "target_3_10001_16000");
        Assert.assertEquals(coverage.targets().size(), bedFeatureCollection.targetCount());
    }

    @Test
    public void testSparkGenomeReadCountsBigBins() throws IOException {
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
        final ReadCountCollection coverage = ReadCountCollectionUtils.parse(outputFile);
        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 4);
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 16000);
        Assert.assertEquals(bedFeatureCollection.target(2).getName(), "target_3_1_16000");
        Assert.assertEquals(coverage.targets().size(), bedFeatureCollection.targetCount());
    }

    @Test
    public void testSparkGenomeReadCountsSmallBins()  throws IOException {
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
        final ReadCountCollection proportionalCoverage = ReadCountCollectionUtils.parse(outputFile);
        Assert.assertTrue(proportionalCoverage.records().stream().anyMatch(t -> Math.abs(t.getDouble(0)) > 1e-10));

        // The reads are all in three bins of contig 3 with values {.5, .25, .25}
        Assert.assertTrue(proportionalCoverage.records().stream().filter(t -> t.getContig().equals("3")).anyMatch(t -> Math.abs(t.getDouble(0)) > .2));
        Assert.assertTrue(Math.abs(proportionalCoverage.records().stream().filter(t -> t.getContig().equals("3")).mapToDouble(t -> t.getDouble(0)).sum() - 1.0) < 1e-10);

        // raw coverage
        final ReadCountCollection coverage = ReadCountCollectionUtils.parse(new File(outputFile.getAbsolutePath() + SparkGenomeReadCounts.RAW_COV_OUTPUT_EXTENSION));
        Assert.assertTrue(coverage.records().stream().anyMatch(t -> Math.abs(t.getDouble(0)) > 1e-10));

        // The reads are all in three bins of contig 3 with values
        Assert.assertEquals(coverage.records().stream().filter(t -> t.getContig().equals("3")).filter(t -> Math.abs(t.getDouble(0)) >= 1).count(), 3);

        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        final TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 16000/2000 * 4); // 4 is the number of contigs in the fasta file
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 4000);
        Assert.assertEquals(bedFeatureCollection.target(2).getName(), "target_1_4001_6000");
        Assert.assertEquals(bedFeatureCollection.target(8).getName(), "target_2_1_2000");
        Assert.assertEquals(bedFeatureCollection.target(17).getName(), "target_3_2001_4000");
        Assert.assertEquals(proportionalCoverage.targets().size(), bedFeatureCollection.targetCount());
    }

    private ReadCountCollection loadReadCountCollection(File outputFile) {
        try {
            return ReadCountCollectionUtils.parse(outputFile);
        } catch (final IOException ioe) {
            Assert.fail("IO Exception in automated test.  Possible misconfiguration?", ioe);
            return null;
        }
    }

    @Test
    public void testSparkGenomeReadCountsInterval() {
        final File outputFile = createTempFile(BAM_FILE.getName(), ".cov");
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BINSIZE_SHORT_NAME, "10000",
                "-L", "1"
        };
        runCommandLine(arguments);

        final ReadCountCollection proportionalCoverage = loadReadCountCollection(outputFile);
        Assert.assertTrue(proportionalCoverage.records().stream().noneMatch(t -> t.getContig().equals("2") || t.getContig().equals("3")));

        // raw coverage
        final ReadCountCollection rawCoverage = loadReadCountCollection(new File(outputFile.getAbsolutePath() + SparkGenomeReadCounts.RAW_COV_OUTPUT_EXTENSION));
        Assert.assertTrue(rawCoverage.records().stream().noneMatch(t -> t.getContig().equals("2") || t.getContig().equals("3")));

        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        final TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertTrue(bedFeatureCollection.targets().stream().allMatch(t -> t.getContig().equals("1")));
    }

    @Test
    public void testSparkGenomeReadCountsTwoIntervals() {
        final File outputFile = createTempFile(BAM_FILE.getName(), ".cov");
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BINSIZE_SHORT_NAME, "10000",
                "-L", "1", "-L", "2"
        };
        runCommandLine(arguments);

        final ReadCountCollection proportionalCoverage = loadReadCountCollection(outputFile);
        Assert.assertTrue(proportionalCoverage.records().stream().anyMatch(t -> t.getContig().equals("1")));
        Assert.assertTrue(proportionalCoverage.records().stream().anyMatch(t -> t.getContig().equals("2")));
        Assert.assertTrue(proportionalCoverage.records().stream().noneMatch(t -> t.getContig().equals("3") || t.getContig().equals("4")));

        // raw coverage
        final ReadCountCollection rawCoverage = loadReadCountCollection(new File(outputFile.getAbsolutePath() + SparkGenomeReadCounts.RAW_COV_OUTPUT_EXTENSION));
        Assert.assertTrue(rawCoverage.records().stream().anyMatch(t -> t.getContig().equals("1")));
        Assert.assertTrue(rawCoverage.records().stream().anyMatch(t -> t.getContig().equals("2")));
        Assert.assertTrue(rawCoverage.records().stream().noneMatch(t -> t.getContig().equals("3") || t.getContig().equals("4")));

        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        final TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertTrue(bedFeatureCollection.targets().stream().allMatch(t -> (t.getContig().equals("1")) || (t.getContig().equals("2"))));
    }

    @Test
    public void testSparkGenomeReadCountsSubContig() {
        final File outputFile = createTempFile(BAM_FILE.getName(), ".cov");
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.OUTPUT_FILE_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BINSIZE_SHORT_NAME, "100",
                "-L", "1:1-500"
        };
        runCommandLine(arguments);

        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        final TargetCollection<BEDFeature> bedFeatureCollection = TargetCollectionUtils.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertTrue(bedFeatureCollection.targets().stream().allMatch(t -> t.getContig().equals("1")));
        Assert.assertEquals(bedFeatureCollection.targets().size(), 5);

        for (int i = 0; i < bedFeatureCollection.targets().size(); i ++) {
            Assert.assertEquals(bedFeatureCollection.targets().get(i).getStart(), i*100 + 1);
            Assert.assertEquals(bedFeatureCollection.targets().get(i).getEnd(), (i+1)*100);
        }
    }
}
