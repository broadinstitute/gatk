package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.pon.coverage.CoveragePoNNormalizationResult;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Stores the results of a tangent normalization.
 */
public final class PCATangentNormalizationResult implements CoveragePoNNormalizationResult {
    /**
     * Name of the column that contains the I.D. of the PoN <i>"eigensample"</i>.
     */
    public static final String PON_SAMPLE_BETA_HAT_COLUMN_NAME = "PON.SAMPLE";

    private final ReadCountCollection tangentNormalized;
    private final ReadCountCollection preTangentNormalized;
    private final RealMatrix tangentBetaHats;
    private final ReadCountCollection targetFactorNormalizedCounts;

    public PCATangentNormalizationResult(final ReadCountCollection tangentNormalized, final ReadCountCollection preTangentNormalized, final RealMatrix tangentBetaHats, final ReadCountCollection targetFactorNormalizedCounts) {
        this.tangentNormalized = tangentNormalized;
        this.preTangentNormalized = preTangentNormalized;
        this.tangentBetaHats = tangentBetaHats;
        this.targetFactorNormalizedCounts = targetFactorNormalizedCounts;
    }

    public ReadCountCollection getTangentNormalized() {
        return tangentNormalized;
    }

    public ReadCountCollection getPreTangentNormalized() {
        return preTangentNormalized;
    }

    public RealMatrix getTangentBetaHats() {
        return tangentBetaHats;
    }

    public ReadCountCollection getTargetFactorNormalizedCounts() {
        return targetFactorNormalizedCounts;
    }

    @Override
    public ReadCountCollection getProfile() {
        return getTangentNormalized();
    }

    /**
     * Write results along with the command line that produced them to specified files.
     * If any file is null, no output will be written for that file.
     */
    public void write(final String commandLine,
                      final File tangentNormalizedFile,
                      final File preTangentNormalizedFile,
                      final File tangentBetaHatsFile,
                      final File targetFactorNormalizedCountsFile) {
        Utils.nonNull(commandLine);
        writeTangentNormalizedOutput(tangentNormalizedFile, commandLine);
        writePreTangentNormalizationOutput(preTangentNormalizedFile, commandLine);
        writeTangentBetaHats(tangentBetaHatsFile, commandLine);
        writeTargetFactorNormalizedOutput(targetFactorNormalizedCountsFile, commandLine);
    }

    /**
     * Writes the tangent-normalized counts in the format described in {@link ReadCountCollectionUtils}.
     */
    private void writeTangentNormalizedOutput(final File file, final String commandLine) {
        writeOutput(file, tangentNormalized, commandLine, "Tangent normalized coverage profile");
    }

    /**
     * Writes the pre-tangent-normalized counts in the format described in {@link ReadCountCollectionUtils}.
     * If the output file is null, we do not write a file.
     */
    private void writePreTangentNormalizationOutput(final File file, final String commandLine) {
        writeOutput(file, preTangentNormalized, commandLine, "Pre tangent normalized coverage profile");
    }

    /**
     * Write the beta-hats.
     * If the output file is null, we do not write a file.
     */
    private void writeTangentBetaHats(final File file, final String commandLine) {
        if (file == null) {
            return;
        }
        final List<String> countColumnNames = targetFactorNormalizedCounts.columnNames();
        final List<String> columnNames = new ArrayList<>(countColumnNames.size() + 1);
        columnNames.add(PON_SAMPLE_BETA_HAT_COLUMN_NAME);
        columnNames.addAll(countColumnNames);
        final TableColumnCollection columns = new TableColumnCollection(columnNames);
        try (final TableWriter<Integer> writer = TableUtils.writer(file, columns,
                (i, dataLine) -> dataLine.append(Integer.toString(i)).append(tangentBetaHats.getRow(i)))) {
            writer.writeComment("fileFormat  = tsv");
            writer.writeComment("commandLine = " + commandLine);
            writer.writeComment("title = Tangent normalization Beta Hats");
            for (int i = 0; i < tangentBetaHats.getRowDimension(); i++) {
                writer.writeRecord(i);
            }
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, ex.getMessage());
        }
    }

    /**
     * Writes the target-factor-normalized counts in the format described in {@link ReadCountCollectionUtils}.
     * If the output file is null, we do not write a file.
     */
    private void writeTargetFactorNormalizedOutput(final File file, final String commandLine) {
        writeOutput(file, targetFactorNormalizedCounts, commandLine, "Target factor normalized target counts");
    }

    private void writeOutput(final File file, final ReadCountCollection counts, final String commandLine, final String title) {
        if (file == null) {
            return;
        }
        try {
            ReadCountCollectionUtils.write(file, counts,
                    "fileFormat = tsv",
                    "commandLine = " + commandLine,
                    "title = " + title);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(file, ex.getMessage());
        }
    }
}
