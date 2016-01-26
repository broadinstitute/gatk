package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.List;
import java.util.SortedMap;
import java.util.stream.Collectors;


public final class TargetCoverageUtils {
    //TODO: these columns disagree with our standard uppercase ones
    //TODO: this is necessary to match gatk-public, but it's really awful and must be fixed
    public static final String TARGET_NAME_COLUMN = "name";
    public static final String CONTIG_COLUMN = "contig";
    public static final String START_COLUMN = "start";
    public static final String END_COLUMN = "stop";

    private TargetCoverageUtils() {}

    /**
     * read a list of targets with coverage from a tab-separated file with header line:
     * name contig  start   stop    <sample name>
     * where the final column is numerical data, generally coverage or normalized coverage.
     *
     * There can only be one sample in the file.
     *
     * NOTE:  For now, the sample name must be the fifth column (index = 4).
     */
    public static List<TargetCoverage> readTargetsWithCoverage(final File targetsFile) {
        try (final TableReader<TargetCoverage> reader = TableUtils.reader(targetsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.containsAll(TARGET_NAME_COLUMN, CONTIG_COLUMN, START_COLUMN, END_COLUMN) || (columns.columnCount() < 5))
                        throw formatExceptionFactory.apply("Bad header");
                    if (columns.columnCount() > 5)
                        throw formatExceptionFactory.apply("Bad header -- more than one sample included in the input file.");

                    //return the lambda to translate dataLines into targets
                    return (dataLine) -> new TargetCoverage(dataLine.get(TARGET_NAME_COLUMN),
                            new SimpleInterval(dataLine.get(CONTIG_COLUMN), dataLine.getInt(START_COLUMN), dataLine.getInt(END_COLUMN)), dataLine.getDouble(4));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(targetsFile, e);
        }
    }

    public static TargetCollection<TargetCoverage> readModeledTargetFileIntoTargetCollection(final File file) {
        try {
            final List<TargetCoverage> targetCoverages = readTargetsWithCoverage(file);
            return new HashedListTargetCollection<>(targetCoverages);
        } catch (final UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(file, e);
        }
    }

    /**
     * write a list of targets with coverage to file, without specifying any comments
     *
     */
    public static void writeTargetsWithCoverage(final File outFile, final String sampleName,
                                                final List<TargetCoverage> targets) {
        writeTargetsWithCoverage(outFile, sampleName, targets, new String [0]);
    }

    /**
     * write a list of targets with coverage to file
     *
     */
    public static void writeTargetsWithCoverage(final File outFile, final String sampleName,
                                                final List<TargetCoverage> targets, final String[] comments) {

        Utils.nonNull(outFile, "Output file cannot be null.");
        Utils.nonNull(sampleName, "Sample name cannot be null.");
        Utils.nonNull(targets, "Targets cannot be null.");
        Utils.nonNull(comments, "Comments cannot be null.");

        final boolean areTargetIntervalsAllPopulated = targets.stream().allMatch(t -> t.getInterval() != null);
        if (!areTargetIntervalsAllPopulated) {
            throw new UserException("Cannot write target coverage file with any null intervals.");
        }

        try (final TableWriter<TargetCoverage> writer = TableUtils.writer(outFile,
                new TableColumnCollection(TARGET_NAME_COLUMN, CONTIG_COLUMN, START_COLUMN, END_COLUMN, sampleName),
                //lambda for filling an initially empty DataLine
                (target, dataLine) -> {
                    final SimpleInterval interval = target.getInterval();
                    final String contig = interval.getContig();
                    final int start = interval.getStart();
                    final int end = interval.getEnd();
                    final String name = target.getName();
                    final double coverage = target.getCoverage();
                    dataLine.append(name, contig).append(start, end).append(coverage);
                })) {

            for (String comment : comments) {
                writer.writeComment(comment);
            }
            for (final TargetCoverage target : targets) {
                writer.writeRecord(target);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }

    /**
     * Write a list of targets with coverage to file with dummy names
     * @param outFile File to write targets with coverage. Never {@code null}
     * @param sampleName Name of sample being written. Never {@code null}
     * @param byKeySorted Map of simple-intervals to their copy-ratio. Never {@code null}
     * @param comments Comments to add to header of coverage file.
     */
    public static <N extends Number> void writeTargetsWithCoverageFromSimpleInterval(final File outFile, final String sampleName,
                                                                  final SortedMap<SimpleInterval, N> byKeySorted,
                                                                  final String[] comments) {

        Utils.nonNull(outFile, "Output file cannot be null.");
        Utils.nonNull(sampleName, "Sample name cannot be null.");
        Utils.nonNull(byKeySorted, "Targets cannot be null.");
        Utils.nonNull(comments, "Comments cannot be null.");

        final boolean areTargetIntervalsAllPopulated = byKeySorted.keySet().stream().allMatch(t -> t != null);
        if (!areTargetIntervalsAllPopulated) {
            throw new UserException("Cannot write target coverage file with any null intervals.");
        }

        try (final TableWriter<TargetCoverage> writer = TableUtils.writer(outFile,
                new TableColumnCollection(TARGET_NAME_COLUMN, CONTIG_COLUMN, START_COLUMN, END_COLUMN, sampleName),
                //lambda for filling an initially empty DataLine
                (target, dataLine) -> {
                    final SimpleInterval interval = target.getInterval();
                    final String contig = interval.getContig();
                    final int start = interval.getStart();
                    final int end = interval.getEnd();
                    final String name = target.getName();
                    final double coverage = target.getCoverage();
                    dataLine.append(name, contig).append(start, end).append(coverage);
                })) {

            for (String comment : comments) {
                writer.writeComment(comment);
            }
            for (final Locatable locatable : byKeySorted.keySet()) {
                final String targetName = createDummyTargetName(locatable);
                final SimpleInterval targetInterval = new SimpleInterval(locatable);
                final double targetCoverage = byKeySorted.get(locatable).doubleValue();
                final TargetCoverage target = new TargetCoverage(targetName, targetInterval, targetCoverage);
                writer.writeRecord(target);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }


    /**
     * Creates a string for a locatable that can be used when creating dummy target names
     * @param locatable The genome region to create a unique dummy target name. Never {@code null}
     * @return never {@code null}
     */
    public static String createDummyTargetName(final Locatable locatable){
        Utils.nonNull(locatable, "Output file cannot be null.");
        return "target_" + locatable.getContig() + "_" + String.valueOf(locatable.getStart()) + "_" + String.valueOf(locatable.getEnd());
    }

    /**
     * write a list of targets without coverage to BED file, only storing chr, start, end, and name
     *
     */
    public static void writeTargetsAsBed(final File outFile, final List<Target> targets) {

        Utils.nonNull(outFile, "Output file cannot be null.");
        Utils.nonNull(targets, "Targets cannot be null.");

        final boolean areTargetIntervalsAllPopulated = targets.stream().allMatch(t -> t.getInterval() != null);
        if (!areTargetIntervalsAllPopulated) {
            throw new UserException("Cannot write target file with any null intervals.");
        }

        try (final FileWriter fw = new FileWriter(outFile)) {
            fw.write("##CONTIG  START\t\tEND\tNAME\n");
            for (Target t: targets) {
                final List<String> fields = Arrays.asList(t.getContig(), String.valueOf(t.getInterval().getGA4GHStart()),
                        String.valueOf(t.getInterval().getGA4GHEnd()), t.getName());

                final String lineString = StringUtils.join(fields, '\t');
                fw.write(lineString + "\n");
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }
}
