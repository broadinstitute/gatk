package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountFileHeaderKey;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.ReadCountType;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Integration test for {@link CollectBinnedReadCounts}
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public class CollectBinnedReadCountsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/collectreadcounts";

    //input files
    private final static File NA12878_BAM = new File(TEST_SUB_DIR, "exome-read-counts-NA12878.bam");
    private final static File NA12778_BAM = new File(TEST_SUB_DIR, "exome-read-counts-NA12778.bam");
    private final static File NA12878_NA12778_BAM_MERGED = new File(TEST_SUB_DIR, "exome-read-counts-NA12878-NA12778-merged.bam");
    private final static File INTERVALS_LIST = new File(TEST_SUB_DIR, "exome-read-counts-intervals.list");
    private final static File INTERVALS_BED = new File(TEST_SUB_DIR, "exome-read-counts-intervals.tsv");
    private final static File INTERVALS_LIST_DUPS = new File(TEST_SUB_DIR, "exome-read-counts-intervals_dups.list");
    private final static File REFERENCE = new File(publicTestDir + "large/human_g1k_v37.20.21.fasta");

    //output files
    private final static File NA12878_SIMPLE_COUNT_EXPECTED_OUTPUT = new File(TEST_SUB_DIR, "");
    private final static File NA12878_GC_BINNING_ONLY_EXPECTED_OUTPUT = new File(TEST_SUB_DIR, "");
    private final static File NA12878_FL_BINNING_ONLY_EXPECTED_OUTPUT = new File(TEST_SUB_DIR, "");
    private final static File NA12878_FL_AND_GC_BINNING_EXPECTED_OUTPUT = new File(TEST_SUB_DIR, "");

    private final static List<String> headerKeysToCompare = Arrays.asList(
            ReadCountFileHeaderKey.READ_COUNT_TYPE.getHeaderKeyName(),
            ReadCountFileHeaderKey.BINNING_CONFIGURATION.getHeaderKeyName());

    @Override
    public String getTestedClassName() {
        return CollectBinnedReadCounts.class.getSimpleName();
    }

    @DataProvider(name="correctRunData")
    public Object[][] correctRunData() {
        return new Object[][] {
//                {
//                        new File[]{ NA12878_BAM },
//                        INTERVALS_LIST,
//                        NA12878_SIMPLE_COUNT_EXPECTED_OUTPUT,
//                        ReadCountType.SIMPLE_COUNT.getReadCountTypeName(),
//                        new String[0]
//                },
//                {
//                        new File[]{ NA12878_BAM },
//                        INTERVALS_LIST,
//                        NA12878_GC_BINNING_ONLY_EXPECTED_OUTPUT,
//                        ReadCountType.BINNED.getReadCountTypeName(),
//                        new String[] { "-" + CollectBinnedReadCounts.INCLUDE_GC_BINS_SHORT_NAME, "true",
//                                "-" + CollectBinnedReadCounts.NUMBER_GC_BINS_SHORT_NAME, "10"}
//                },
//                {
//                        new File[]{ NA12878_BAM },
//                        INTERVALS_LIST,
//                        NA12878_FL_BINNING_ONLY_EXPECTED_OUTPUT,
//                        ReadCountType.BINNED.getReadCountTypeName(),
//                        new String[] { "-" + CollectBinnedReadCounts.INCLUDE_FRAGMENT_LENGTH_BINS_SHORT_NAME, "true",
//                                "-" + CollectBinnedReadCounts.FRAGMENT_LENGTH_MIN_BIN_VALUE_SHORT_NAME, "0",
//                                "-" + CollectBinnedReadCounts.FRAGMENT_LENGTH_MAX_BIN_VALUE_SHORT_NAME, "1000",
//                                "-" + CollectBinnedReadCounts.NUMBER_FRAGMENT_LENGTH_BINS_SHORT_NAME, "20"}
//                },
                {
                        new File[]{ NA12878_BAM },
                        INTERVALS_LIST,
                        NA12878_FL_AND_GC_BINNING_EXPECTED_OUTPUT,
                        REFERENCE,
                        ReadCountType.BINNED.getReadCountTypeName(),
                        new String[] { "-" + CollectBinnedReadCounts.INCLUDE_GC_BINS_SHORT_NAME, "true",
                                "-" + CollectBinnedReadCounts.NUMBER_GC_BINS_SHORT_NAME, "100",
                                "-" + CollectBinnedReadCounts.INCLUDE_FRAGMENT_LENGTH_BINS_SHORT_NAME, "true",
                                "-" + CollectBinnedReadCounts.FRAGMENT_LENGTH_MIN_BIN_VALUE_SHORT_NAME, "0",
                                "-" + CollectBinnedReadCounts.FRAGMENT_LENGTH_MAX_BIN_VALUE_SHORT_NAME, "500",
                                "-" + CollectBinnedReadCounts.NUMBER_FRAGMENT_LENGTH_BINS_SHORT_NAME, "1",
                                "-" + IntervalArgumentCollection.INTERVAL_MERGING_RULE_LONG_NAME, IntervalMergingRule.OVERLAPPING_ONLY.toString()}
                }

        };
    }

    //TODO add test that checks that exception is thrown when no interval file is supplied

    //TODO add test that checks that exception is thrown when a target file is supplied instead of picard style intervals

    //TODO add test that checks for validity of covariate configuration parameters

    //TODO test that inteval merging rule is set correctly
    @Test(dataProvider = "correctRunData")
    public void testCorrectRun(final File[] bamFiles, final File intervalFile, final File expectedOutputFile,
                               final File reference, final String readCountType, final String[] additionalArguments) {
        final File outputFile = BaseTest.createTempFile("GATK4-collectCounts-outputTable", ".tmp");
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        ));
        if (intervalFile != null) {
            arguments.add("-L");
            arguments.add(intervalFile.getAbsolutePath());
        }
        if (reference != null) {
            arguments.add("-R");
            arguments.add(reference.getAbsolutePath());
        }
        arguments.add("-" + CollectBinnedReadCounts.READ_COUNT_TYPE_SHORT_NAME);
        arguments.add(readCountType);
        Arrays.asList(bamFiles).forEach(bam -> {
            arguments.add("-" + StandardArgumentDefinitions.INPUT_SHORT_NAME);
            arguments.add(bam.getAbsolutePath());
        });
        arguments.addAll(Arrays.asList(additionalArguments));
        System.err.println("COMMAND LINE: " + Arrays.toString(arguments.toArray()));
        runCommandLine(arguments);
        IntegerReadCountFileComparator.AssertEquals(outputFile, expectedOutputFile, headerKeysToCompare);
    }

    private static class IntegerReadCountFileComparator {
        private final List<int[]> values;
        private final int columnCount;
        private final int rowCount;
        private final String[] columnNames;
        private final Map<String, String> headerValuesMap;

        private final static String TAB_SEPARATOR_REGEX = "\\t";
        private final static String COMMENT_PREFIX_START_REGEX = "^#.*$";

        /**
         * Private constructor
         *
         * @param file input file
         * @throws IOException if there is any problems while reading the file
         * @throws UserException.BadInput if non-integer value is encountered below column names line or file
         * is not formatted correctly
         */
        private IntegerReadCountFileComparator(final File file) throws IOException {
            Utils.nonNull(file, "the file cannot be null");
            final BufferedReader lineReader = new BufferedReader(new FileReader(file));
            headerValuesMap = new HashMap<>(10);
            String lastLine;

            //read in the header values
            String headerKeyValuePattern = null;
            while ((lastLine = lineReader.readLine()) != null) {
                if (lastLine.matches(COMMENT_PREFIX_START_REGEX)) {
                    headerValuesMap.put(lastLine.replaceAll(headerKeyValuePattern, "$1"), lastLine.replaceAll(headerKeyValuePattern, "$2"));
                } else {
                    break;
                }
            }

            //read in table's column values
            columnNames = lastLine == null ? null : lineReader.readLine().split(TAB_SEPARATOR_REGEX);
            Utils.nonNull(columnNames, "the table does not have a header");
            columnCount = columnNames.length;

            //read in table's integer values, throw exception if encounter a non-integer type
            values = new ArrayList<>();
            int currentRowCount = 0;
            while ((lastLine = lineReader.readLine()) != null) {
                final String[] lineValues = lastLine.split(TAB_SEPARATOR_REGEX);
                if (lineValues.length != columnCount) {
                    throw new UserException.BadInput(String.format("line %d of the table has different number of values than the column name line", currentRowCount + 1));
                }
                values.add(new int[columnCount]);
                for (int index = 0; index < columnCount; index ++) {
                    try {
                        values.get(currentRowCount)[index] = Integer.parseInt(lineValues[index]);
                    } catch (NumberFormatException ex) {
                        throw new UserException.BadInput(String.format("Table contains value %s which is non-integer ", lineValues[index]));
                    }
                }
                currentRowCount++;
            }

            rowCount = currentRowCount;
            lineReader.close();
        }

        private static void AssertEquals(final File f1, final File f2, final List<String> headerKeysToCompare) {
            try {
                IntegerReadCountFileComparator firstTable = new IntegerReadCountFileComparator(f1);
                IntegerReadCountFileComparator secondTable = new IntegerReadCountFileComparator(f2);
                //check that mandatory header key value pairs are the same.
                headerKeysToCompare.stream().forEach(key -> {
                    Assert.assertTrue(firstTable.headerValuesMap.containsKey(key));
                    Assert.assertTrue(secondTable.headerValuesMap.containsKey(key));
                    Assert.assertEquals(firstTable.headerValuesMap.get(key), secondTable.headerValuesMap.get(key));
                });
                //check that tables have the same dimensions and columns
                Assert.assertEquals(firstTable.columnCount,secondTable.columnCount, "Different number of columns");
                Assert.assertEquals(firstTable.rowCount,secondTable.rowCount, "Different number of rows");
                Assert.assertEquals(firstTable.columnNames,secondTable.columnNames, "Different column names");
                //check that tables have same values row by row
                IntStream.range(0, firstTable.rowCount).forEach(row -> Assert.assertEquals(firstTable.values.get(row), secondTable.values.get(row)));
            } catch (final IOException ex) {
                Assert.fail("Failed to read actual or expected table files ", ex);
            }
        }
    }

}