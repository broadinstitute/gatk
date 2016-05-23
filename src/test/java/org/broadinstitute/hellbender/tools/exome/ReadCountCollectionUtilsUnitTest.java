package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.bed.BEDCodec;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.Collectors;

/**
 * Unit tests for {@link ReadCountCollectionUtils}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadCountCollectionUtilsUnitTest {

    private static final File TEST_FILE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/exome");
    private static final File FULL_CORRECT_FILE = new File(TEST_FILE_DIR, "rcc-test-full-counts.txt");

    private static final String CONTIG_START_END = TargetTableColumns.CONTIG.toString()
            + "\t" + TargetTableColumns.START.toString() + "\t" + TargetTableColumns.END.toString();

    private static final String CONTIG_START_END_NAME = CONTIG_START_END + "\t" + TargetTableColumns.NAME.toString();

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullFile() throws IOException {
        ReadCountCollectionUtils.parse(null);
    }

    @Test
    public void testReadFullCorrectFile() throws IOException {
        ReadCountCollectionUtils.parse(FULL_CORRECT_FILE);
    }

    @Test(expectedExceptions = IOException.class)
    public void testReadNonExistingFile() throws IOException {
        ReadCountCollectionUtils.parse(new File(TEST_FILE_DIR, "false-" + FULL_CORRECT_FILE.getName()));
    }

    @Test(expectedExceptions = IOException.class)
    public void testReadFromDirectory() throws IOException {
        ReadCountCollectionUtils.parse(TEST_FILE_DIR);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testEmptyFile() throws IOException {
        final File emptyFile = createTempFile();
        ReadCountCollectionUtils.parse(emptyFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testCommentHeaderOnly() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoCountColumns() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME);
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testNoSpecialColumns() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println("SAMPLE1\tSAMPLE2");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRepeatedSampleNames() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE1");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRepeatedSpecialColumnNames() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println("CONTIG\tCONTIG\tEND\tNAME\tSAMPLE1\tSAMPLE1");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadEmptyCounts() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE2");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testRepeatedTargets() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE2");
        writer.println("1\t100\t200\ttgt_0\t1.1\t2.2");
        writer.println("1\t201\t300\ttgt_1\t1.1\t2.2");
        writer.println("2\t400\t500\ttgt_0\t-1.1E-7\t-2.2E-8");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test
    public void testReadFullFormattedFile() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println(CONTIG_START_END_NAME + "\tSAMPLE1\tSAMPLE2");
        writer.println("1\t100\t200\ttgt_0\t1.1\t2.2");
        writer.println("2\t200\t300\ttgt_1\t-1.1E-7\t-2.2E-8");
        writer.close();
        final ReadCountCollection subject = ReadCountCollectionUtils.parse(testFile);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.columnNames(), Arrays.asList("SAMPLE1", "SAMPLE2"));
        Assert.assertEquals(subject.targets().stream().map(Target::getName).collect(Collectors.toList()), Arrays.asList("tgt_0", "tgt_1"));
        Assert.assertEquals(subject.targets().stream().map(Target::getInterval).collect(Collectors.toList()), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("2", 200, 300)));
        Assert.assertEquals(subject.targets().size(), 2);
        final RealMatrix counts = subject.counts();
        Assert.assertEquals(counts.getRowDimension(), 2);
        Assert.assertEquals(counts.getColumnDimension(), 2);
        Assert.assertEquals(counts.getEntry(0, 0), 1.1, 0.0001);
        Assert.assertEquals(counts.getEntry(0, 1), 2.2, 0.0001);
        Assert.assertEquals(counts.getEntry(1, 0), -1.1E-7, 0.000000001);
        Assert.assertEquals(counts.getEntry(1, 1), -2.2E-8, 0.000000001);
    }

    @Test
    public void testReadTargetNameOnlyFormattedFile() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("## comment 1");
        writer.println("## comment 2");
        writer.println("SAMPLE2\tSAMPLE1\t" + TargetTableColumns.NAME.toString());
        writer.println("1.1\t2.2\ttgt_0");
        writer.println("-1.1E-7\t-2.2E-8\ttgt_1");
        writer.close();
        final ReadCountCollection subject = ReadCountCollectionUtils.parse(testFile);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.columnNames(), Arrays.asList("SAMPLE2", "SAMPLE1"));
        Assert.assertEquals(subject.targets().stream().map(Target::getName).collect(Collectors.toList()), Arrays.asList("tgt_0", "tgt_1"));
        Assert.assertEquals(subject.targets().stream().map(Target::getInterval).collect(Collectors.toList()), Arrays.asList(null, null));
        Assert.assertEquals(subject.targets().size(), 2);
        final RealMatrix counts = subject.counts();
        Assert.assertEquals(counts.getRowDimension(), 2);
        Assert.assertEquals(counts.getColumnDimension(), 2);
        Assert.assertEquals(counts.getEntry(0, 0), 1.1, 0.0001);
        Assert.assertEquals(counts.getEntry(0, 1), 2.2, 0.0001);
        Assert.assertEquals(counts.getEntry(1, 0), -1.1E-7, 0.000000001);
        Assert.assertEquals(counts.getEntry(1, 1), -2.2E-8, 0.000000001);
    }

    @Test
    public void testReadIntervalsOnlyFile() throws IOException {
        final File targetFile = createTempFile();
        final PrintWriter targetWriter = new PrintWriter(targetFile);
        targetWriter.println("#" + CONTIG_START_END_NAME);
        targetWriter.println("1\t99\t200\tTGT_0");
        targetWriter.println("2\t199\t300\tTGT_1");
        targetWriter.close();
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("contig\tSAMPLE3\tstop\tSAMPLE2\tstart");
        writer.println("1\t1.1\t200\t2.2\t100");
        writer.println("2\t-1.1E-7\t300\t-2.2E-8\t200");
        writer.close();

        final ReadCountCollection subject = ReadCountCollectionUtils.parse(testFile,
                TargetCollectionUtils.fromBEDFeatureFile(targetFile, new BEDCodec()), false);
        Assert.assertNotNull(subject);
        Assert.assertEquals(subject.columnNames(), Arrays.asList("SAMPLE3", "SAMPLE2"));
        Assert.assertEquals(subject.targets().stream().map(Target::getInterval).collect(Collectors.toList()), Arrays.asList(new SimpleInterval("1", 100, 200), new SimpleInterval("2", 200, 300)));
        Assert.assertEquals(subject.targets().stream().map(Target::getName).collect(Collectors.toList()), Arrays.asList("TGT_0", "TGT_1"));
        Assert.assertEquals(subject.targets().size(), 2);
        final RealMatrix counts = subject.counts();
        Assert.assertEquals(counts.getRowDimension(), 2);
        Assert.assertEquals(counts.getColumnDimension(), 2);
        Assert.assertEquals(counts.getEntry(0, 0), 1.1, 0.0001);
        Assert.assertEquals(counts.getEntry(0, 1), 2.2, 0.0001);
        Assert.assertEquals(counts.getEntry(1, 0), -1.1E-7, 0.000000001);
        Assert.assertEquals(counts.getEntry(1, 1), -2.2E-8, 0.000000001);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testOneTooFewValueInLine() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("contig\tSAMPLE3\tstop\tSAMPLE2\tstart");
        writer.println("1.1\t200\t2.2\t100");
        writer.println("2\t-1.1E-7\t300\t-2.2E-8\t200");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testOneTooManyValueInLine() throws IOException {
        final File testFile = createTempFile();
        final PrintWriter writer = new PrintWriter(testFile);
        writer.println("contig\tSAMPLE3\tstop\tSAMPLE2\tstart");
        writer.println("1\t1.1\t200\t2.2\t100\t21.0");
        writer.println("2\t-1.1E-7\t300\t-2.2E-8\t200");
        writer.close();
        ReadCountCollectionUtils.parse(testFile);
    }

    private File createTempFile() throws IOException {
        final File result = File.createTempFile("file", ".test");
        result.deleteOnExit();
        return result;
    }
}
