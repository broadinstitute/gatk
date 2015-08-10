package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
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
 * Created by David Benjamin 7/15/15
 */
public final class SegmentUtils {
    /**
     * read a list of intervals without calls from a segfile
     */
    public static List<SimpleInterval> readIntervalsFromSegfile(final File segmentsFile) {
        try (final TableReader<SimpleInterval> reader = TableUtils.reader(segmentsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesAll(0, "Sample", "Chromosome", "Start", "End")) {//ignore last two columns: NTARGETS and AVERAGE
                        throw formatExceptionFactory.apply("Bad header");
                    }

                    // return the lambda to translate dataLines into uncalled segments.
                    return (dataLine) ->
                            new SimpleInterval(dataLine.get(1), dataLine.getInt(2), dataLine.getInt(3));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(segmentsFile, e);
        }
    }

    /**
     * read a list of intervals with calls from a file
     */
    public static List<CalledInterval> readCalledIntervalsFromSegfile(final File segmentsFile) {
        try (final TableReader<CalledInterval> reader = TableUtils.reader(segmentsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("Sample", "Chromosome", "Start", "End", "Call")) {
                        throw formatExceptionFactory.apply("Bad header");
                    }
                    // return the lambda to translate dataLines into called segments.
                    return (dataLine) -> new CalledInterval(
                            new SimpleInterval(dataLine.get(1), dataLine.getInt(2), dataLine.getInt(3)), dataLine.get(4));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(segmentsFile, e);
        }
    }

    /**
     * write a list of intervals with calls to file
     */
    public static void writeCalledIntervalsToSegfile(final File outFile, List<CalledInterval> segments, final String sample) {
        try (final TableWriter<CalledInterval> writer = TableUtils.writer(outFile,
                    new TableColumnCollection("Sample", "Chromosome", "Start", "End", "Call"),

                    //lambda for filling an initially empty DataLine
                    (ci, dataLine) -> {
                        dataLine.append(sample, ci.getContig()).append(ci.getStart(), ci.getEnd())
                                .append(ci.getCall());
                    })) {
            for (final CalledInterval ci : segments) {
                if (ci == null) {
                    throw new IllegalArgumentException("Segments list contains a null.");
                }
                writer.writeRecord(ci);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }

    /**
     * the mean of all overlapping targets' coverages
     *
     * @throws IllegalStateException if overlapping targets have not been assigned or if no overlapping targets were found.
     */
    public static double meanTargetCoverage(final Locatable loc, final TargetCollection<TargetCoverage> targets) {
        Utils.nonNull(loc, "Can't get mean coverage of null segment.");
        Utils.nonNull(targets, "Mean target coverage requires non-null targets collection.");

        final List<TargetCoverage> myTargets = targets.targets(loc);

        if (myTargets.size() == 0) {
            throw new IllegalStateException("Empty segment -- no overlapping targets.");
        }
        return myTargets.stream().mapToDouble(TargetCoverage::getCoverage).average().getAsDouble();
    }
}
