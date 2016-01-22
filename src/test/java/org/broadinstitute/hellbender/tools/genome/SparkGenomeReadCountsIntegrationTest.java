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
        List<TargetCoverage> targetCoverages = TargetCoverageUtils.readTargetsWithCoverage(outputFile);
        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        TargetCollection<BEDFeature> bedFeatureCollection = TargetCollections.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 8);
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 16000);
        Assert.assertEquals(bedFeatureCollection.target(5).getName(), "target_3_10001_16000");
        Assert.assertEquals(targetCoverages.size(), bedFeatureCollection.targetCount());
    }

    @Test
    public void testSparkGenomeReadCountsBigBins() {
        final File outputFile = createTempFile(BAM_FILE.getName(),".cov");
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
        List<TargetCoverage> targetCoverages = TargetCoverageUtils.readTargetsWithCoverage(outputFile);
        final File bedFile = new File(outputFile.getAbsolutePath()+".bed");
        Assert.assertTrue(bedFile.exists());
        Assert.assertTrue(bedFile.length() > 0);

        TargetCollection<BEDFeature> bedFeatureCollection = TargetCollections.fromBEDFeatureFile(bedFile, new BEDCodec());
        Assert.assertEquals(bedFeatureCollection.targets().size(), 4);
        Assert.assertEquals(bedFeatureCollection.target(1).getEnd(), 16000);
        Assert.assertEquals(bedFeatureCollection.target(2).getName(), "target_3_1_16000");
        Assert.assertEquals(targetCoverages.size(), bedFeatureCollection.targetCount());
    }

    @Test
    public void testDropContigsFromSequence() throws IOException{
        final Set<String> testContigsToDrop = new HashSet<>(Arrays.asList("3"));
        ReferenceFileSource testFasta = new ReferenceFileSource(REFERENCE_FILE.getAbsolutePath());
        final SAMSequenceDictionary originalSequenceDictionary = testFasta.getReferenceSequenceDictionary(null);
        SAMSequenceDictionary modifiedSequenceDictionary = originalSequenceDictionary;
        modifiedSequenceDictionary = SparkGenomeReadCounts.dropContigsFromSequence(originalSequenceDictionary, testContigsToDrop);
        SAMSequenceDictionary finalSequenceDictionary = modifiedSequenceDictionary;
        Assert.assertEquals(finalSequenceDictionary.size(), 3);

    }
}
