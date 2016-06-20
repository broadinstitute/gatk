package org.broadinstitute.hellbender.tools.exome;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.List;

/**
 * Target table file reader.
 *
 * <p>
 * <a href="#format"/>
 * The minimum content of a target table file consists of 4 columns indicate the target's name, that must be unique across
 * targets and its coordinates: enclosing contig, start and end positions. The start and end positions are inclusive
 * and 1-based (the first position index is 1).
 * </p>
 * <p>
 * Example:
 * <pre>
 *         CONTIG   START   END     NAME
 *         1        100     200     target_0
 *         1        300     400     target_1
 *         2        100     500     target_2
 * </pre>
 * </p>
 * <p>
 * There might be additional columns indicating arbitrary target annotations such as per-sample coverage, gc content
 * and so forth.
 * </p>
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetTableReader extends TableReader<Target> {
    protected final static Logger logger = LogManager.getLogger(TargetTableReader.class);
    private TargetTableAnnotationManager annotationCollection;

    public TargetTableReader(final Reader reader) throws IOException {
        super(reader);
    }

    public TargetTableReader(final File file) throws IOException {
        super(file);
    }

    /**
     * Read targets from a target file in the format of {@link TargetWriter}.
     *
     * @param targetsFile Target table file.
     * @return never {@code null}
     */
    public static List<Target> readTargetFile(final File targetsFile) {
        Utils.regularReadableUserFile(targetsFile);
        logger.log(Level.INFO, String.format("Reading targets from file '%s' ...", targetsFile.getAbsolutePath()));
        final List<Target> inputTargets;
        try (final TargetTableReader reader = new TargetTableReader(targetsFile)) {
            inputTargets = reader.toList();
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read input targets file", ex);
        }
        return inputTargets;
    }

    @Override
    protected void processColumns(final TableColumnCollection columns) {
        TableUtils.checkMandatoryColumns(columns, TargetTableColumn.MANDATORY_COLUMNS, this::formatException);
        annotationCollection = new TargetTableAnnotationManager(getSource(), columns);
    }

    @Override
    protected Target createRecord(final DataLine dataLine) {
        final String contig = dataLine.get(TargetTableColumn.CONTIG);
        final int start = dataLine.getInt(TargetTableColumn.START);
        final int end = dataLine.getInt(TargetTableColumn.END);
        if (start < 0) {
            throw formatException("the start position must not be negative: " + start);
        } else if (start > end) {
            throw formatException(String.format("the start position %d cannot be greater than the end position %d", start, end));
        }
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        final String name = dataLine.get(TargetTableColumn.NAME);
        return new Target(name, interval, annotationCollection.createTargetAnnotationCollection(dataLine));
    }
}