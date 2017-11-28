package org.broadinstitute.hellbender.tools.genome;

import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetTableReader;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class SparkGenomeReadCountsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/genome");
    private static final File BAM_FILE = new File(TEST_FILE_DIR, "HCC1143_chr3_1K_11K.tiny.bam");
    private static final File REFERENCE_FILE = new File("src/test/resources/hg19mini.fasta");
    private static final String TSV_EXT = SparkGenomeReadCounts.TSV_EXT;

    @Test
    public void testSparkGenomeReadCounts() throws IOException {
        final File outputFile = createTempFile(BAM_FILE.getName(), TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "10000",
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        final ReadCountCollection rcc = ReadCountCollectionUtils.parse(outputFile);

        final File intervalsFile = new File(FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + SparkGenomeReadCounts.INTERVALS_EXT);
        Assert.assertTrue(intervalsFile.exists());
        Assert.assertTrue(intervalsFile.length() > 0);
        final List<Target> intervals = TargetTableReader.readTargetFile(intervalsFile);
        Assert.assertEquals(intervals.size(), 8);
        Assert.assertEquals(intervals.get(1).getEnd(), 16000);
        Assert.assertEquals(intervals.get(5).getName(), "target_3_10001_16000");
        Assert.assertEquals(rcc.targets().size(), intervals.size());
    }

    @Test
    public void testSparkGenomeReadCountsBigBins() throws IOException {
        final File outputFile = createTempFile(BAM_FILE.getName(), TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "16000",
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        Assert.assertTrue(outputFile.length() > 0);
        final ReadCountCollection rcc = ReadCountCollectionUtils.parse(outputFile);

        final File intervalsFile = new File(FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + SparkGenomeReadCounts.INTERVALS_EXT);
        Assert.assertTrue(intervalsFile.exists());
        Assert.assertTrue(intervalsFile.length() > 0);
        final List<Target> intervals = TargetTableReader.readTargetFile(intervalsFile);
        Assert.assertEquals(intervals.size(), 4);
        Assert.assertEquals(intervals.get(1).getEnd(), 16000);
        Assert.assertEquals(intervals.get(2).getName(), "target_3_1_16000");
        Assert.assertEquals(rcc.targets().size(), intervals.size());
    }

    @Test
    public void testSparkGenomeReadCountsSmallBins()  throws IOException {
        final File outputFile = createTempFile(BAM_FILE.getName(), TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "2000",
                "-" + SparkGenomeReadCounts.WRITE_HDF5_SHORT_NAME, "false"
        };
        runCommandLine(arguments);
        Assert.assertTrue(outputFile.exists());
        final File hdf5File = new File(FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + SparkGenomeReadCounts.HDF5_EXT);
        Assert.assertFalse(hdf5File.exists());
        Assert.assertTrue(outputFile.length() > 0);

        final ReadCountCollection rcc = ReadCountCollectionUtils.parse(outputFile);
        Assert.assertTrue(rcc.records().stream().anyMatch(t -> Math.abs(t.getDouble(0)) > 1e-10));

        // The reads are all in three bins of contig 3 with values
        Assert.assertEquals(rcc.records().stream().filter(t -> t.getContig().equals("3")).filter(t -> Math.abs(t.getDouble(0)) >= 1).count(), 3);

        final File intervalsFile = new File(FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + SparkGenomeReadCounts.INTERVALS_EXT);
        Assert.assertTrue(intervalsFile.exists());
        Assert.assertTrue(intervalsFile.length() > 0);
        final List<Target> intervals = TargetTableReader.readTargetFile(intervalsFile);
        Assert.assertEquals(intervals.size(), 16000/2000 * 4); // 4 is the number of contigs in the fasta file
        Assert.assertEquals(intervals.get(1).getEnd(), 4000);
        Assert.assertEquals(intervals.get(2).getName(), "target_1_4001_6000");
        Assert.assertEquals(intervals.get(8).getName(), "target_2_1_2000");
        Assert.assertEquals(intervals.get(17).getName(), "target_3_2001_4000");
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
        final File outputFile = createTempFile(BAM_FILE.getName(), TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "10000",
                "-L", "1"
        };
        runCommandLine(arguments);

        final ReadCountCollection rcc = loadReadCountCollection(outputFile);
        Assert.assertTrue(rcc.records().stream().noneMatch(t -> t.getContig().equals("2") || t.getContig().equals("3")));

        final File intervalsFile = new File(FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + SparkGenomeReadCounts.INTERVALS_EXT);
        Assert.assertTrue(intervalsFile.exists());
        Assert.assertTrue(intervalsFile.length() > 0);
        final List<Target> intervals = TargetTableReader.readTargetFile(intervalsFile);
        Assert.assertTrue(intervals.stream().allMatch(t -> t.getContig().equals("1")));
    }

    @Test
    public void testSparkGenomeReadCountsTwoIntervals() {
        final File outputFile = createTempFile(BAM_FILE.getName(), TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "10000",
                "-L", "1", "-L", "2"
        };
        runCommandLine(arguments);

        final ReadCountCollection rcc = loadReadCountCollection(outputFile);
        Assert.assertTrue(rcc.records().stream().anyMatch(t -> t.getContig().equals("1")));
        Assert.assertTrue(rcc.records().stream().anyMatch(t -> t.getContig().equals("2")));
        Assert.assertTrue(rcc.records().stream().noneMatch(t -> t.getContig().equals("3") || t.getContig().equals("4")));

        final File intervalsFile = new File(FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + SparkGenomeReadCounts.INTERVALS_EXT);
        Assert.assertTrue(intervalsFile.exists());
        Assert.assertTrue(intervalsFile.length() > 0);
        final List<Target> intervals = TargetTableReader.readTargetFile(intervalsFile);
        Assert.assertTrue(intervals.stream().allMatch(t -> (t.getContig().equals("1")) || (t.getContig().equals("2"))));
    }

    @Test
    public void testSparkGenomeReadCountsSubContig() {
        final File outputFile = createTempFile(BAM_FILE.getName(), TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "100",
                "-L", "1:1-500"
        };
        runCommandLine(arguments);

        final File intervalsFile = new File(FilenameUtils.removeExtension(outputFile.getAbsolutePath()) + SparkGenomeReadCounts.INTERVALS_EXT);
        Assert.assertTrue(intervalsFile.exists());
        Assert.assertTrue(intervalsFile.length() > 0);
        final List<Target> intervals = TargetTableReader.readTargetFile(intervalsFile);
        Assert.assertTrue(intervals.stream().allMatch(t -> t.getContig().equals("1")));
        Assert.assertEquals(intervals.size(), 5);

        for (int i = 0; i < intervals.size(); i ++) {
            Assert.assertEquals(intervals.get(i).getStart(), i*100 + 1);
            Assert.assertEquals(intervals.get(i).getEnd(), (i+1)*100);
        }
    }

    @Test
    public void testSparkGenomeReadCountsHdf5Writing() throws IOException {
        final File tsvFile = createTempFile(BAM_FILE.getName(),TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tsvFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "10000",
                "-" + SparkGenomeReadCounts.WRITE_HDF5_SHORT_NAME
        };
        runCommandLine(arguments);
        Assert.assertTrue(tsvFile.exists());
        final File hdf5File = new File(FilenameUtils.removeExtension(tsvFile.getAbsolutePath()) + SparkGenomeReadCounts.HDF5_EXT);
        Assert.assertTrue(hdf5File.exists());
        Assert.assertTrue(hdf5File.length() > 0);

        final SimpleCountCollection scc = SimpleCountCollection.read(hdf5File);
        final ReadCountCollection rccTsv = ReadCountCollectionUtils.parse(tsvFile);

        Assert.assertEquals(scc.getCounts(), rccTsv.counts().transpose().getRow(0));
        Assert.assertEquals(scc.getSampleName(), rccTsv.columnNames().get(0));
        Assert.assertEquals(scc.getIntervals(), rccTsv.targets().stream().map(Target::getInterval).collect(Collectors.toList()));

        // Make sure we are putting integer counts in the HDF5
        Assert.assertEquals(MathUtils.sum(scc.getCounts()), 4.0);
    }

    @Test
    public void testSparkGenomeReadCountsProportionalCoverageWriting() throws IOException {
        final File tsvFile = createTempFile(BAM_FILE.getName(),TSV_EXT);
        final String[] arguments = {
                "--disableSequenceDictionaryValidation",
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REFERENCE_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, BAM_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, tsvFile.getAbsolutePath(),
                "-" + SparkGenomeReadCounts.BIN_LENGTH_SHORT_NAME, "10000",
                "-" + SparkGenomeReadCounts.WRITE_PROPORTIONAL_COVERAGE_SHORT_NAME
        };
        runCommandLine(arguments);
        Assert.assertTrue(tsvFile.exists());
        final File pCovFile = new File(FilenameUtils.removeExtension(tsvFile.getAbsolutePath()) + SparkGenomeReadCounts.PROPORTIONAL_COVERAGE_EXT);
        Assert.assertTrue(pCovFile.exists());
        Assert.assertTrue(pCovFile.length() > 0);

        final ReadCountCollection rccPCov = ReadCountCollectionUtils.parse(pCovFile);
        final ReadCountCollection rccTsv = ReadCountCollectionUtils.parse(tsvFile);

        final double totalCounts = MathUtils.sum(rccTsv.counts().getColumn(0));
        final double[] pCovCalculatedFromTsv = Arrays.stream(rccTsv.counts().getColumn(0)).map(c -> c / totalCounts).toArray();
        Assert.assertEquals(rccPCov.counts().getColumn(0), pCovCalculatedFromTsv);
        Assert.assertEquals(rccPCov.columnNames().get(0), rccTsv.columnNames().get(0));
        Assert.assertEquals(rccPCov.targets(), rccTsv.targets());

        // Make sure we are putting proportional coverage in the pcov file
        Assert.assertEquals(MathUtils.sum(rccPCov.counts().getColumn(0)), 1.0);
    }
}