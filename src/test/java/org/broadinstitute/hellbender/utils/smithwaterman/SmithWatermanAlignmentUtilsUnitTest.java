package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.stream.Collectors;

public final class SmithWatermanAlignmentUtilsUnitTest extends GATKBaseTest {
    private static final String SMITH_WATERMAN_PARAMETERS_TABLE_COLUMN_HEADER_LINE =
            Arrays.stream(SmithWatermanAlignmentUtils.SmithWatermanParametersTableReader.SmithWatermanParametersTableColumn.values())
                    .map(SmithWatermanAlignmentUtils.SmithWatermanParametersTableReader.SmithWatermanParametersTableColumn::toString)
                    .collect(Collectors.joining("\t")) + '\n';
    private static final SWParameters EXPECTED_SMITH_WATERMAN_PARAMETERS = new SWParameters(100, -200, -300, -10);

    private static final File VALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE = createTempFile("smith-waterman-parameters-table", ".tsv");
    private static final File INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITH_MORE_THAN_ONE_ROW_OF_VALUES = createTempFile("invalid-smith-waterman-parameters-table", ".tsv");
    private static final File INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITHOUT_COLUMN_HEADER_LINE = createTempFile("invalid-smith-waterman-parameters-table", ".tsv");

    @BeforeTest
    private static void setup() {
        try {
            Files.write(
                    VALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE.toPath(),
                    SMITH_WATERMAN_PARAMETERS_TABLE_COLUMN_HEADER_LINE.getBytes());
            Files.write(
                    VALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE.toPath(),
                    String.format("%d\t%d\t%d\t%d\n",
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMatchValue(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMismatchPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapOpenPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapExtendPenalty()).getBytes(),
                    StandardOpenOption.APPEND);

            Files.write(
                    INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITH_MORE_THAN_ONE_ROW_OF_VALUES.toPath(),
                    SMITH_WATERMAN_PARAMETERS_TABLE_COLUMN_HEADER_LINE.getBytes());
            Files.write(
                    INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITH_MORE_THAN_ONE_ROW_OF_VALUES.toPath(),
                    String.format("%d\t%d\t%d\t%d\n",
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMatchValue(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMismatchPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapOpenPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapExtendPenalty()).getBytes(),
                    StandardOpenOption.APPEND);
            Files.write(
                    INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITH_MORE_THAN_ONE_ROW_OF_VALUES.toPath(),
                    String.format("%d\t%d\t%d\t%d\n",
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMatchValue(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMismatchPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapOpenPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapExtendPenalty()).getBytes(),
                    StandardOpenOption.APPEND);

            Files.write(
                    INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITHOUT_COLUMN_HEADER_LINE.toPath(),
                    String.format("%d\t%d\t%d\t%d\n",
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMatchValue(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getMismatchPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapOpenPenalty(),
                            EXPECTED_SMITH_WATERMAN_PARAMETERS.getGapExtendPenalty()).getBytes(),
                    StandardOpenOption.APPEND);
        } catch (final IOException e) {
            throw new GATKException.ShouldNeverReachHereException("Unable to write test files.");
        }
    }

    @Test
    public void testReadSmithWatermanParametersFromTSV() {
        final SWParameters swParameters = SmithWatermanAlignmentUtils.readSmithWatermanParametersFromTSV(
                new GATKPath(VALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE.toURI().toString()));
        assertSWParametersEqual(swParameters, EXPECTED_SMITH_WATERMAN_PARAMETERS);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadSmithWatermanParametersFromTSVWithMoreThanOneRowOfValues() {
        SmithWatermanAlignmentUtils.readSmithWatermanParametersFromTSV(
                new GATKPath(INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITH_MORE_THAN_ONE_ROW_OF_VALUES.toURI().toString()));
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadSmithWatermanParametersFromTSVWithoutColumnHeaderLine() {
        SmithWatermanAlignmentUtils.readSmithWatermanParametersFromTSV(
                new GATKPath(INVALID_SMITH_WATERMAN_PARAMETERS_TABLE_FILE_WITHOUT_COLUMN_HEADER_LINE.toURI().toString()));
    }

    private static void assertSWParametersEqual(final SWParameters result,
                                                final SWParameters expected) {
        Assert.assertEquals(result.getMatchValue(), expected.getMatchValue());
        Assert.assertEquals(result.getMismatchPenalty(), expected.getMismatchPenalty());
        Assert.assertEquals(result.getGapOpenPenalty(), expected.getGapOpenPenalty());
        Assert.assertEquals(result.getGapExtendPenalty(), expected.getGapExtendPenalty());
    }
}