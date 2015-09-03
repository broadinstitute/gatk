package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Tests for FilterReads
 */
public final class FilterReadsIntegrationTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/FilterReads/");

    @Override
    public String getTestedClassName() {
        return FilterReads.class.getSimpleName();
    }

    @DataProvider(name="filterTestsData")
    public Object[][] filterTestData() {
        return new Object[][]{
                // inputFileName, outputExtension, filter, readListFile, sortOrder, writeReadsFile, expectedReadsCount, expectedOutputCount
                {"all_aligned.sam", ".sam", "includeAligned", null, "unsorted", true, 4, 4},
                {"all_aligned.sam", ".sam", "excludeAligned", null, "unsorted", false, 0, 0},

                {"mixed_aligned.sam", ".sam", "includeAligned", null, "unsorted", false, 0, 8},
                {"mixed_aligned.sam", ".sam", "excludeAligned", null, "unsorted", false, 0, 2},
                {"mixed_aligned.sam", ".sam", "includeAligned", null, "coordinate", true, 10, 8},
                {"mixed_aligned.sam", ".sam", "excludeAligned", null, "queryname", false, 0, 2},
                {"mixed_aligned.sam", ".sam", "includeReadList", "readlist.txt", "unsorted", true, 10, 2},
                {"mixed_aligned.sam", ".sam", "excludeReadList", "readlist.txt", "unsorted", false, 0, 8},

                {"unmapped.sam", ".sam", "includeAligned", null, "unsorted", false, 0, 0},
                {"unmapped.sam", ".sam", "excludeAligned", null, "unsorted", false, 0, 10},
                {"unmapped.sam", ".sam", "includeAligned", null, "coordinate", false, 0, 0},
                {"unmapped.sam", ".sam", "excludeAligned", null, "queryname", true, 10, 10},
                {"unmapped.sam", ".sam", "includeReadList", "readlist.txt", "coordinate", false, 0, 2},
                {"unmapped.sam", ".sam", "excludeReadList", "readlist.txt", "queryname", false, 0, 8},

                {"all_aligned.bam", ".bam", "includeAligned", null, "unsorted", true, 4, 4},
                {"all_aligned.bam", ".bam", "excludeAligned", null, "unsorted", false, 0, 0},

                {"mixed_aligned.bam", ".bam", "includeAligned", null, "unsorted", true, 10, 8},
                {"mixed_aligned.bam", ".bam", "excludeAligned", null, "unsorted", false, 0, 2},
                {"mixed_aligned.bam", ".bam", "includeAligned", null, "coordinate", false, 0, 8},
                {"mixed_aligned.bam", ".bam", "excludeAligned", null, "queryname", false, 0, 2},
                {"mixed_aligned.bam", ".bam", "includeReadList", "readlist.txt", "unsorted", false, 0, 2},
                {"mixed_aligned.bam", ".bam", "excludeReadList", "readlist.txt", "unsorted", false, 0, 8},

                {"unmapped.bam", ".bam", "includeAligned", null, "unsorted", false, 0, 0},
                {"unmapped.bam", ".bam", "excludeAligned", null, "unsorted", false, 0, 10},
                {"unmapped.bam", ".bam", "includeAligned", null, "coordinate", true, 10, 0},
                {"unmapped.bam", ".bam", "excludeAligned", null, "queryname", false, 0, 10},
                {"unmapped.bam", ".bam", "includeReadList", "readlist.txt", "coordinate", false, 0, 2},
                {"unmapped.bam", ".bam", "excludeReadList", "readlist.txt", "queryname", false, 0, 8},
        };
    }

    @Test(dataProvider="filterTestsData")
    public void testReadFilter(
            final String inputFileName,
            final String outputExtension,
            final String filter,
            final String readListFile,
            final String sortOrder,
            final boolean writeReadsFile,
            final int expectedReadsCount,
            final int expectedOutputCount) throws Exception {
        final List<String> args = new ArrayList<>();
        final File inputFile = new File(TEST_DATA_DIR, inputFileName);

        args.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
        args.add(inputFile.getAbsolutePath());

        final String outFileName = BaseTest.createTempFile(inputFileName, outputExtension).getAbsolutePath();
        args.add("-"+ StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        args.add(outFileName);

        if (null != filter) {
            args.add("--FILTER");
            args.add(filter);
        }
        if (null != readListFile) {
            final String readlistFile = new File(TEST_DATA_DIR, readListFile).getAbsolutePath();
            args.add("--READ_LIST_FILE");
            args.add(readlistFile);
        }

        if (!writeReadsFile) { // defaults to true
            args.add("--WRITE_READS_FILE");
            args.add("false");
        }

        switch (sortOrder) {
            case "coordinate": {
                args.add("--SORT_ORDER");
                args.add("coordinate");
                break;
            }
            case "queryname": {
                args.add("--SORT_ORDER");
                args.add("queryname");
                break;
            }
            case "unsorted":
            default:
                break;
        }

        Assert.assertNull(runCommandLine(args));

        Assert.assertEquals(getReadCounts(outFileName), expectedOutputCount);
        Assert.assertTrue(validateSortOrder(outFileName, sortOrder));
        Assert.assertEquals(getReadsFileCount(inputFile, new File(outFileName), writeReadsFile), expectedReadsCount);
    }

    private int getReadCounts(final String resultFileName) {
        final File path = new File(resultFileName);
        IOUtil.assertFileIsReadable(path);
        final SamReader in = SamReaderFactory.makeDefault().open(path);
        int count = 0;
        for (@SuppressWarnings("unused") final SAMRecord rec : in) {
            count++;
        }
        CloserUtil.close(in);
        return count;
    }

    private boolean validateSortOrder(final String resultFileName, final String sortOrderName) throws IOException {
        final File path = new File(resultFileName);
        IOUtil.assertFileIsReadable(path);

        try (final SamReader in = SamReaderFactory.makeDefault().open(path);) {
            final SAMFileHeader header = in.getFileHeader();
            final SAMFileHeader.SortOrder hdrOrder = header.getSortOrder();
            switch (sortOrderName) {
                case "coordinate": {
                    return hdrOrder == SAMFileHeader.SortOrder.coordinate;
                }
                case "queryname": {
                    return hdrOrder == SAMFileHeader.SortOrder.queryname;
                }
                case "unsorted":
                default: // don't care what it is
                    return true;
            }
        }
    }

    private int getReadsFileCount(final File inputFile, final File outputFile, final boolean writeReadsFile) throws IOException {
        int count = 0;
        if (writeReadsFile) {
            final File readsFile = new File(outputFile.getParentFile(), IOUtil.basename(inputFile) + ".reads");
            try (final InputStream inStream = IOUtil.openFileForReading(readsFile);
                  final BufferedLineReader br = new BufferedLineReader(inStream)) {
                while (null != br.readLine()) {
                    count++;
                }
            }
        }
        return count;
    }

}
