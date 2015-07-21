package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Interval;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Simple data structure to pass and read/write List of AllelicCounts.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountCollection {
    private final List<AllelicCount> counts;

    public AllelicCountCollection() {
        this.counts = new ArrayList<>();
    }

    /**
     * Constructor that reads (sequence, position, reference count, alternate count) from the specified file.
     * @param inputFile     file to read from
     */
    public AllelicCountCollection(final File inputFile) {
        this.counts = new ArrayList<>();
        try (final TableReader<AllelicCount> reader = TableUtils.reader(inputFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("SEQ", "POS", "REF_COUNT", "ALT_COUNT"))
                        throw formatExceptionFactory.apply("Bad header");

                    // return the lambda to translate dataLines into AllelicCounts.
                    return (dataLine) -> {
                        final Interval interval = new Interval(dataLine.get(0), dataLine.getInt(1), dataLine.getInt(1));
                        final int refReadCount = dataLine.getInt(2);
                        final int altReadCount = dataLine.getInt(3);
                        return new AllelicCount(interval, refReadCount, altReadCount);
                    };
                })) {
            for (final AllelicCount count : reader) {
                this.counts.add(count);
            }
        } catch (final Exception e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e.getMessage());
        }
    }

    /**
     * Adds (interval, reference count, alternate count) to respective lists.
     * @param interval      site in 1-based format
     * @param refReadCount  number of reads at site matching the reference
     * @param altReadCount  number of reads at site different from the reference
     */
    public void add(final Interval interval, final int refReadCount, final int altReadCount) {
        final AllelicCount count = new AllelicCount(interval, refReadCount, altReadCount);
        counts.add(count);
    }

    /** Returns a list of the allelic counts.   */
    public List<AllelicCount> getCounts() {
        return Collections.unmodifiableList(counts);
    }

    /**
     * Writes out (sequence, position, reference count, alternate count) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection("SEQ", "POS", "REF_COUNT", "ALT_COUNT"),
                //lambda for filling an initially empty DataLine
                (count, dataLine) -> {
                    final Interval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    dataLine.append(interval.getContig()).append(interval.getEnd(), refReadCount, altReadCount);
                })) {
            for (final AllelicCount count : counts) {
                writer.writeRecord(count);
            }
        } catch (final Exception e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e.getMessage());
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCountCollection)) {
            return false;
        }

        final AllelicCountCollection allelicCountCollection = (AllelicCountCollection) o;
        return counts.equals(allelicCountCollection.counts);
    }

    @Override
    public int hashCode() {
        return counts.hashCode();
    }
}
