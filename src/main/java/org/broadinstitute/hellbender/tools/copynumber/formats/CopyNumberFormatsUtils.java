package org.broadinstitute.hellbender.tools.copynumber.formats;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public final class CopyNumberFormatsUtils {
    public static final String COMMENT_PREFIX = "@";    //SAMTextHeaderCodec.HEADER_LINE_START; we need TableReader to treat SAM header as comment lines
    public static final String DOUBLE_FORMAT = "%.6f";

    private CopyNumberFormatsUtils() {}

    public static String formatDouble(final double value) {
        return String.format(DOUBLE_FORMAT, value);
    }

    /**
     * Extracts column names from a TSV file
     */
    public static TableColumnCollection readColumnsFromHeader(final File inputFile) {
        IOUtils.canReadFile(inputFile);
        List<String> columns = null;
        try (final XReadLines reader = new XReadLines(inputFile)) {
            while (reader.hasNext()) {
                String nextLine = reader.next();
                if (!nextLine.startsWith(COMMENT_PREFIX)) {
                    columns = Arrays.asList(nextLine.split(TableUtils.COLUMN_SEPARATOR_STRING));
                    break;
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile);
        }
        if (columns == null) {
            throw new UserException.BadInput(String.format(
                    "The input file %s does not have a header (starting with comment character %s).",
                    inputFile.getAbsolutePath(), COMMENT_PREFIX));
        }
        if (columns.stream().distinct().count() != columns.size()) {
            throw new UserException.BadInput("Column headers must all be unique.");
        }
        return new TableColumnCollection(columns);
    }
}
