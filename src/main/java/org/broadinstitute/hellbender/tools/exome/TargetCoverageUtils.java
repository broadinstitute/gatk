package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by davidben on 7/17/15.
 */
public final class TargetCoverageUtils {
    /**
     * read a list of targets with coverage from a file
     */
    public static List<TargetCoverage> readTargetsWithCoverage(final File file) throws IOException {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(file,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesAll(0, "CONTIG", "START", "END", "NAME"))   //coverage is fifth column w/ header = <sample name>
                        throw formatExceptionFactory.apply("Bad header");

                    // return the lambda to translate dataLines into targets.
                    return (dataLine) -> new TargetCoverage(dataLine.get(3),
                            new SimpleInterval(dataLine.get(0), dataLine.getInt(1), dataLine.getInt(2)),
                            dataLine.getDouble(4));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (UncheckedIOException e) {
            throw e.getCause();
        }
    }
}
