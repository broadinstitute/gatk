package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Simple data structure to pass and write pulldown results.  Should probably replace with a more generic class later.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Pulldown {
    private int size;
    private final IntervalList intervals;
    private final List<Integer> refReadCounts;
    private final List<Integer> altReadCounts;

    public Pulldown(final SAMFileHeader header) {
        Utils.nonNull(header, "SAMFileHeader must be supplied.");

        this.intervals = new IntervalList(header);
        this.refReadCounts = new ArrayList<>();
        this.altReadCounts = new ArrayList<>();
    }

    /**
     * Constructor that reads (sequence, position, reference count, alternate count) from the specified file and
     * uses external SAMFile header to construct Pulldown.
     * @param inputFile     file to read from
     * @param header        SAMFile header for IntervalList
     * TODO remove dependency from IntervalList/SamLocusIterator on external header once LocusWalker implemented?
     */
    public Pulldown(final File inputFile, final SAMFileHeader header) {
        this(header);

        try (final TableReader<List<String>> reader = TableUtils.reader(inputFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("SEQ", "POS", "REF", "ALT"))
                        throw formatExceptionFactory.apply("Bad header");

                    // return the lambda to translate dataLines into Pulldown rows.
                    return (dataLine) -> Arrays.asList(dataLine.get(0), dataLine.get(1),
                            dataLine.get(2), dataLine.get(3));
                })) {
            for (final List<String> pulldownRow : reader) {
                final String seq = pulldownRow.get(0);
                final int pos = Integer.parseInt(pulldownRow.get(1));
                final int refReadCount = Integer.parseInt(pulldownRow.get(2));
                final int altReadCount = Integer.parseInt(pulldownRow.get(3));
                this.add(new Interval(seq, pos, pos), refReadCount, altReadCount);
            }
        } catch (final Exception e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e.getMessage());
        }
    }

    /**
     * Adds (interval, reference count, alternate count) to respective lists.
     * @param interval      heterozygous SNP site in the format (seq, pos, pos)
     * @param refReadCount  number of reads at SNP site matching the reference
     * @param altReadCount  number of reads at SNP site different from the reference
     */
    public void add(final Interval interval, final int refReadCount, final int altReadCount) {
        intervals.add(interval);
        refReadCounts.add(refReadCount);
        altReadCounts.add(altReadCount);
        size++;
    }

    /** Returns the IntervalList of heterozygous SNP sites.   */
    public IntervalList getIntervals() { return IntervalList.copyOf(intervals);     }

    /**
     * Writes out (sequence, position, reference count, alternate count) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile) {

        try (final TableWriter<Integer> writer = TableUtils.writer(outputFile,
                new TableColumnCollection("SEQ", "POS", "REF", "ALT"),
                //lambda for filling an initially empty DataLine
                (index, dataLine) -> {
                    final Interval interval = intervals.getIntervals().get(index);
                    final int refReadCount = refReadCounts.get(index);
                    final int altReadCount = altReadCounts.get(index);
                    dataLine.append(interval.getContig()).append(interval.getEnd())
                            .append(refReadCount).append(altReadCount);
                })) {
            for (int index = 0; index < size; index++) {
                writer.writeRecord(index);
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
        if (!(o instanceof Pulldown)) {
            return false;
        }

        final Pulldown pulldown = (Pulldown) o;
        return intervals.equals(pulldown.intervals) && refReadCounts.equals(pulldown.refReadCounts)
                && altReadCounts.equals(pulldown.altReadCounts);
    }
}
