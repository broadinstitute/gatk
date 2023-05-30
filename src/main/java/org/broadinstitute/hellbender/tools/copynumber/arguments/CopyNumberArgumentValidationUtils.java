package org.broadinstitute.hellbender.tools.copynumber.arguments;

import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryUtils;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.DetermineGermlineContigPloidy;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractLocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AbstractRecordCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AnnotatedInterval;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.OptionalInt;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyNumberArgumentValidationUtils {
    private static final Logger logger = LogManager.getLogger(CopyNumberArgumentValidationUtils.class);

    private CopyNumberArgumentValidationUtils() {}

    /**
     * Validate that the interval-argument collection parameters minimally modify the input intervals.
     */
    public static void validateIntervalArgumentCollection(final IntervalArgumentCollection intervalArgumentCollection) {
        Utils.validateArg(intervalArgumentCollection.getIntervalSetRule() == IntervalSetRule.UNION,
                "Interval set rule must be set to UNION.");
        Utils.validateArg(intervalArgumentCollection.getIntervalExclusionPadding() == 0,
                "Interval exclusion padding must be set to 0.");
        Utils.validateArg(intervalArgumentCollection.getIntervalPadding() == 0,
                "Interval padding must be set to 0.");
        Utils.validateArg(intervalArgumentCollection.getIntervalMergingRule() == IntervalMergingRule.OVERLAPPING_ONLY,
                "Interval merging rule must be set to OVERLAPPING_ONLY.");
    }

    /**
     * Validate that a list of locatables is valid and sorted according to a sequence dictionary and contains no duplicates or overlaps.
     */
    public static <T extends Locatable> void validateIntervals(final List<T> intervals,
                                                               final SAMSequenceDictionary sequenceDictionary) {
        Utils.nonNull(intervals);
        Utils.nonNull(sequenceDictionary);
        Utils.validateArg(intervals.stream().allMatch(i -> IntervalUtils.intervalIsOnDictionaryContig(new SimpleInterval(i), sequenceDictionary)),
                "Records contained at least one interval that did not validate against the sequence dictionary.");
        if (!Ordering.from(IntervalUtils.getDictionaryOrderComparator(sequenceDictionary)).isStrictlyOrdered(intervals)) {
            throw new IllegalArgumentException("Records were not strictly sorted in dictionary order.");
        }
        final OptionalInt failureIndex = IntStream.range(1, intervals.size())
                .filter(i -> IntervalUtils.overlaps(intervals.get(i - 1), intervals.get(i)))
                .findFirst();
        if (failureIndex.isPresent()) {
            final int index = failureIndex.getAsInt();
            throw new IllegalArgumentException(
                    String.format("Records contain at least two overlapping intervals: %s and %s",
                            intervals.get(index - 1), intervals.get(index)));
        }
    }

    /**
     * Compares two non-null sequence dictionaries using sequence index, name, and length only.
     * Less stringent than {@link SAMSequenceDictionary#isSameDictionary}.
     * Use {@link #getValidatedSequenceDictionary} to validate multiple sequence dictionaries from {@link LocatableMetadata}.
     */
    public static boolean isSameDictionary(final SAMSequenceDictionary dictionary1,
                                           final SAMSequenceDictionary dictionary2) {
        Utils.nonNull(dictionary1);
        Utils.nonNull(dictionary2);
        if (dictionary1 == dictionary2) {
            return true;
        }
        final boolean checkContigOrdering = true;
        return SAMSequenceDictionaryUtils.compareDictionaries(dictionary1, dictionary2, checkContigOrdering) ==
                SAMSequenceDictionaryUtils.SequenceDictionaryCompatibility.IDENTICAL;
    }

    /**
     * Resolve intervals from an {@link IntervalArgumentCollection} and a read-count path.
     * If intervals are not specified in the {@link IntervalArgumentCollection}, they are taken from the
     * read-count path.  The sequence dictionary is taken from the read-count path.  A {@link SimpleIntervalCollection}
     * constructed using these intervals and sequence dictionary is returned and can be used for further validation.
     */
    public static SimpleIntervalCollection resolveIntervals(final String readCountPath,
                                                            final IntervalArgumentCollection intervalArgumentCollection,
                                                            final Logger logger) {
        IOUtils.assertFileIsReadable(IOUtils.getPath(readCountPath));
        Utils.nonNull(intervalArgumentCollection);
        Utils.nonNull(logger);

        if (intervalArgumentCollection.intervalsSpecified()) {
            logger.info("Intervals specified...");
            validateIntervalArgumentCollection(intervalArgumentCollection);
        } else {
            logger.info(String.format("Retrieving intervals from read-count file (%s)...", readCountPath));
        }

        final SimpleCountCollection readCounts = BucketUtils.isGcsUrl(readCountPath)
                ? SimpleCountCollection.readFromGCS(readCountPath)
                : SimpleCountCollection.read(new File(readCountPath));
        final SAMSequenceDictionary sequenceDictionary = readCounts.getMetadata().getSequenceDictionary();
        final LocatableMetadata metadata = new SimpleLocatableMetadata(sequenceDictionary);
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified()
                ? intervalArgumentCollection.getIntervals(sequenceDictionary)
                : readCounts.getIntervals();

        return new SimpleIntervalCollection(metadata, intervals);
    }

    /**
     * For all non-null inputs, validate that all metadata are identical and return the metadata.
     */
    @SafeVarargs
    @SuppressWarnings({"varargs"})
    public static <METADATA extends Metadata> METADATA getValidatedMetadata(final AbstractRecordCollection<METADATA, ?> ... recordCollections) {
        Utils.nonNull(recordCollections);
        final Set<METADATA> metadataSet = Stream.of(recordCollections)
                .filter(Objects::nonNull)
                .map(AbstractRecordCollection::getMetadata)
                .collect(Collectors.toSet());
        Utils.nonEmpty(metadataSet, "At least one collection must be non-null.");
        Utils.validateArg(metadataSet.size() == 1, "Metadata do not match.");
        return metadataSet.stream().findFirst().get();
    }

    /**
     * For all non-null inputs, validate that all sequence dictionaries match (using {@link #isSameDictionary})
     * and return the sequence dictionary; otherwise, emit a warning.
     */
    public static SAMSequenceDictionary getValidatedSequenceDictionary(final AbstractLocatableCollection<?, ?> ... locatableCollections) {
        Utils.nonNull(locatableCollections);
        final List<SAMSequenceDictionary> sequenceDictionaries = Stream.of(locatableCollections)
                .filter(Objects::nonNull)
                .map(AbstractLocatableCollection::getMetadata)
                .map(LocatableMetadata::getSequenceDictionary)
                .collect(Collectors.toList());
        Utils.nonEmpty(sequenceDictionaries, "At least one collection must be non-null.");
        if (!IntStream.range(0, sequenceDictionaries.size() - 1).
                allMatch(i -> CopyNumberArgumentValidationUtils.isSameDictionary(sequenceDictionaries.get(i), sequenceDictionaries.get(i + 1)))) {
            logger.warn("Sequence dictionaries do not match across all inputs.");
        }
        return sequenceDictionaries.get(0);
    }

    /**
     * Common method for subsetting and validating read counts in both {@link DetermineGermlineContigPloidy}
     * and {@link GermlineCNVCaller}.
     * @param inputReadCountPaths   for indexed read counts given by GCS paths, counts will be streamed
     * @param specifiedIntervals    intervals to query and subset
     */
    public static Stream<SimpleCountCollection> streamOfSubsettedAndValidatedReadCounts(final List<String> inputReadCountPaths,
                                                                                        final SimpleIntervalCollection specifiedIntervals,
                                                                                        final Logger logger) {
        Utils.nonEmpty(inputReadCountPaths);
        Utils.nonNull(specifiedIntervals);
        Utils.nonNull(logger);
        final int numSamples = inputReadCountPaths.size();
        final Set<SimpleInterval> intervalSubset = new HashSet<>(specifiedIntervals.getRecords());                       //for subsetting local files
        final List<SimpleInterval> mergedIntervalSubset = IntervalUtils.getIntervalsWithFlanks(
                specifiedIntervals.getRecords(), 0, specifiedIntervals.getMetadata().getSequenceDictionary());  //for subsetting GCS files

        return IntStream.range(0, inputReadCountPaths.size()).boxed()
                .map(sampleIndex -> {
                    final String inputReadCountPath = inputReadCountPaths.get(sampleIndex);
                    logger.info(String.format("Aggregating read-count file %s (%d / %d)",
                            inputReadCountPath, sampleIndex + 1, numSamples));
                    final SimpleCountCollection subsetReadCounts = BucketUtils.isGcsUrl(inputReadCountPath)
                            ? SimpleCountCollection.readOverlappingSubsetFromGCS(inputReadCountPath, mergedIntervalSubset)
                            : SimpleCountCollection.readAndSubset(new File(inputReadCountPath), intervalSubset);
                    if (!CopyNumberArgumentValidationUtils.isSameDictionary(
                            subsetReadCounts.getMetadata().getSequenceDictionary(),
                            specifiedIntervals.getMetadata().getSequenceDictionary())) {
                        logger.warn("Sequence dictionary for read-count file {} does not match that " +
                                "in other read-count files.", inputReadCountPath);
                    }
                    Utils.validateArg(subsetReadCounts.size() == intervalSubset.size(),
                            String.format("Intervals for read-count file %s do not contain all specified intervals.",
                                    inputReadCountPath));
                    return subsetReadCounts;
                });
    }

    /**
     * Checks equality of the sequence dictionary and intervals contained in an {@code locatableCollection}
     * against those contained in an {@link AnnotatedIntervalCollection} represented by {@code annotatedIntervalsFile}.
     * If the latter is {@code null}, then {@code null} is returned; otherwise,
     * the {@link AnnotatedIntervalCollection} represented by {@code inputAnnotatedIntervalsFile} is returned
     * if the intervals are equal, and an exception is thrown if they are not.
     */
    public static AnnotatedIntervalCollection validateAnnotatedIntervals(final File annotatedIntervalsFile,
                                                                         final AbstractLocatableCollection<?, ?> locatableCollection,
                                                                         final Logger logger) {
        Utils.nonNull(locatableCollection);
        Utils.nonNull(logger);
        if (annotatedIntervalsFile == null) {
            logger.info("No annotated intervals were provided...");
            return null;
        }
        logger.info("Reading and validating annotated intervals...");
        final AnnotatedIntervalCollection annotatedIntervals = new AnnotatedIntervalCollection(annotatedIntervalsFile);
        final SAMSequenceDictionary sequenceDictionary = locatableCollection.getMetadata().getSequenceDictionary();
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(annotatedIntervals.getMetadata().getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in annotated-intervals file does not match the master sequence dictionary.");
        }
        Utils.validateArg(annotatedIntervals.getIntervals().equals(locatableCollection.getIntervals()),
                "Annotated intervals do not match provided intervals.");
        return annotatedIntervals;
    }

    /**
     * Same as {@link #validateAnnotatedIntervals}, except we only require that {@code annotatedIntervalsFile}
     * contains as a subset all the intervals contained in {@code locatableCollection} along with equality of the sequence dictionaries.
     * The corresponding subset of annotated intervals is returned if appropriate.
     */
    public static AnnotatedIntervalCollection validateAnnotatedIntervalsSubset(final File annotatedIntervalsFile,
                                                                               final AbstractLocatableCollection<?, ?> locatableCollection,
                                                                               final Logger logger) {
        Utils.nonNull(locatableCollection);
        Utils.nonNull(logger);
        if (annotatedIntervalsFile == null) {
            logger.info("No annotated intervals were provided...");
            return null;
        }
        logger.info("Reading and validating annotated intervals...");
        IOUtils.canReadFile(annotatedIntervalsFile);
        final AnnotatedIntervalCollection annotatedIntervals = new AnnotatedIntervalCollection(annotatedIntervalsFile);
        final SAMSequenceDictionary sequenceDictionary = locatableCollection.getMetadata().getSequenceDictionary();
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(annotatedIntervals.getMetadata().getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in annotated-intervals file does not match the master sequence dictionary.");
        }
        final Set<SimpleInterval> intervalsSubset = new HashSet<>(locatableCollection.getIntervals());
        final List<AnnotatedInterval> subsetAnnotatedIntervals = annotatedIntervals.getRecords().stream()
                .filter(i -> intervalsSubset.contains(i.getInterval()))
                .collect(Collectors.toList());
        Utils.validateArg(subsetAnnotatedIntervals.size() == intervalsSubset.size(),
                "Annotated intervals do not contain all specified intervals.");
        return new AnnotatedIntervalCollection(locatableCollection.getMetadata(), subsetAnnotatedIntervals);
    }

    /**
     * Validate that input files and/or directories are readable if they are not {@code null} (i.e., optional inputs).
     */
    public static void validateInputs(final File ... inputs) {
        if (inputs != null) {
            for (final File input : inputs) {
                if (input != null) {
                    if (input.isFile()) {
                        IOUtils.canReadFile(input);
                    } else if (input.isDirectory() && !input.canRead()) {
                        throw new UserException.CouldNotReadInputFile(input.getAbsolutePath());
                    }
                }
            }
        }
    }

    /**
     * Validate that input paths are readable if they are not {@code null} (i.e., optional inputs).
     */
    public static void validateInputs(final String ... inputs) {
        if (inputs != null) {
            for (final String input : inputs) {
                if (input != null) {
                    IOUtils.assertFileIsReadable(IOUtils.getPath(input));
                }
            }
        }
    }

    /**
     * Validate that output files are writeable, whether or not they already exist.
     */
    public static void validateOutputFiles(final File ... outputFiles) {
        Utils.nonNull(outputFiles);
        for (final File outputFile : outputFiles) {
            Utils.nonNull(outputFile);
            if ((outputFile.exists() && !outputFile.canWrite()) || (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException.CouldNotCreateOutputFile(outputFile.getAbsolutePath(), ": The output file is not writeable.");
            }
        }
    }

    /**
     * Validate that output directories are writeable.  If a directory does not exist, create it.
     */
    public static void validateAndPrepareOutputDirectories(final File ... outputDirectories) {
        Utils.nonNull(outputDirectories);
        for (final File outputDirectory : outputDirectories) {
            Utils.nonNull(outputDirectory);
            if (outputDirectory.exists()) {
                if (!outputDirectory.canWrite()) {
                    throw new UserException.CouldNotCreateOutputFile(outputDirectory.getAbsolutePath(), ": The output directory is not writeable.");
                }
            } else {
                try {
                    IOUtils.createDirectory(outputDirectory.getAbsolutePath());
                } catch (final IOException e) {
                    throw new UserException.CouldNotCreateOutputFile(outputDirectory.getAbsolutePath(), ": The output directory does not exist and could not be created.");
                }
            }
        }
    }

    /**
     * File paths that are passed to {@link PythonScriptExecutor} must be canonical (rather than absolute).
     * See https://github.com/broadinstitute/gatk/issues/4724.
     */
    public static String getCanonicalPath(final File file) {
        Utils.nonNull(file);
        try {
            return file.getCanonicalPath();
        } catch (final IOException e) {
            throw new UserException.BadInput(String.format("Could not resolve a canonical file path: %s", file));
        }
    }

    /**
     * File paths that are passed to {@link PythonScriptExecutor} must be canonical (rather than absolute).
     * See https://github.com/broadinstitute/gatk/issues/4724.
     */
    public static String getCanonicalPath(final String filename) {
        Utils.nonEmpty(filename);
        return getCanonicalPath(new File(filename));
    }

    public static String addTrailingSlashIfNecessary(final String outputDir) {
        Utils.nonEmpty(outputDir);
        return outputDir.endsWith(File.separator) ? outputDir : outputDir + File.separator;
    }
}
