package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by davidben on 7/17/15.
 */
public final class TargetCoverageUtils {
    private static final String CONTIG_COLUMN_NAME = "CONTIG";
    private static final String START_COLUMN_NAME = "START";
    private static final String END_COLUMN_NAME = "END";
    private static final String TARGET_NAME_COLUMN_NAME = "NAME";

    /**
     * read a list of targets with coverage from a file
     */
    public static List<TargetCoverage> readTargetsWithCoverage(final File targetsFile) {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(targetsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesAll(0,
                            CONTIG_COLUMN_NAME, START_COLUMN_NAME, END_COLUMN_NAME, TARGET_NAME_COLUMN_NAME))
                            //coverage is fifth column w/ header = <sample name>
                        throw formatExceptionFactory.apply("Bad header");
                    //return the lambda to translate dataLines into targets
                    return (dataLine) -> new TargetCoverage(dataLine.get(3),
                            new SimpleInterval(dataLine.get(0), dataLine.getInt(1), dataLine.getInt(2)),
                            dataLine.getDouble(4));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(targetsFile, e);
        }
    }

    /**
     * read a list of targets with coverage from a file into a TargetCollection
     */
    public static TargetCollection<TargetCoverage> readTargetsWithCoverageIntoTargetCollection(final File file) throws IOException {
        final List<TargetCoverage> targetList = readTargetsWithCoverage(file);
        return new HashedListTargetCollection<>(targetList);
    }

    public static TargetCollection<TargetCoverage> readModeledTargetFileIntoTargetCollection(final File file) throws IOException {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(file,
                (columns, formatExceptionFactory) -> {
                    if (!columns.containsAll("name", "contig", "start", "stop") || (columns.columnCount() < 5))

                        throw formatExceptionFactory.apply("Bad header");
                    //return the lambda to translate dataLines into targets
                    //coverage is fifth column w/ header = <sample name>, so we use the column index.
                    return (dataLine) -> new TargetCoverage(dataLine.get("name"),
                            new SimpleInterval(dataLine.get("contig"), dataLine.getInt("start"), dataLine.getInt("stop")),
                            dataLine.getDouble(4));
                })) {
            return new HashedListTargetCollection<>(reader.stream().collect(Collectors.toList()));
        } catch (final UncheckedIOException e) {
            throw e.getCause();
        }
    }


    /**
     * write a list of targets with coverage to file
     */
    public static void writeTargetsWithCoverage(final File outFile, final String sampleName,
                                                final List<TargetCoverage> targets) {
        try (final TableWriter<TargetCoverage> writer = TableUtils.writer(outFile,
                new TableColumnCollection(
                        CONTIG_COLUMN_NAME, START_COLUMN_NAME, END_COLUMN_NAME, TARGET_NAME_COLUMN_NAME, sampleName),
                //lambda for filling an initially empty DataLine
                (target, dataLine) -> {
                    final SimpleInterval interval = target.getInterval();
                    final String contig = interval.getContig();
                    final int start = interval.getStart();
                    final int end = interval.getEnd();
                    final String name = target.getName();
                    final double coverage = target.getCoverage();
                    dataLine.append(contig).append(start, end).append(name).append(coverage);
                })) {
            for (final TargetCoverage target : targets) {
                writer.writeRecord(target);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }
}
