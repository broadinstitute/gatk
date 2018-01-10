package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.Lazy;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.List;

/**
 * Helper class for {@link SimpleCountCollection} used to read/write HDF5.
 * Class is only visible so that it can be referenced in documentation.
 *
 * <p>
 *     Data is stored in the following HDF5 paths:
 * </p>
 * <ul>
 *     <li>
 *         sample name: /sample_metadata/sample_name
 *     </li>
 *     <li>
 *         sequence dictionary: /locatable_metadata/sequence_dictionary
 *     </li>
 *     <li>
 *         intervals: /intervals
 *     </li>
 *     <li>
 *         counts: /counts/values
 *     </li>
 * </ul>
 * <p>
 *     See {@link HDF5Utils#writeIntervals} for details on the representation of intervals.
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class HDF5SimpleCountCollection {
    private static final String SAMPLE_NAME_PATH = "/sample_metadata/sample_name";
    private static final String SEQUENCE_DICTIONARY_PATH = "/locatable_metadata/sequence_dictionary";
    private static final String INTERVALS_GROUP_NAME = "/intervals";
    private static final String COUNTS_PATH = "/counts/values";

    private final HDF5File file;
    private final Lazy<String> sampleName;
    private final Lazy<SAMSequenceDictionary> sequenceDictionary;
    private final Lazy<List<SimpleInterval>> intervals;
    private final Lazy<RealMatrix> counts;

    /**
     * DEV NOTE: If you are adding attributes that are neither RealMatrix nor a primitive,
     * you must follow the pattern in the constructor (i.e. the Lazy loading pattern).
     * Otherwise, some operations will hang.
     */
    HDF5SimpleCountCollection(final HDF5File file) {
        Utils.nonNull(file, "The input file cannot be null.");
        this.file = file;
        sampleName = new Lazy<>(() -> file.readStringArray(SAMPLE_NAME_PATH)[0]);
        sequenceDictionary = new Lazy<>(() -> {
            final String sequenceDictionaryString = file.readStringArray(SEQUENCE_DICTIONARY_PATH)[0];
            return new SAMTextHeaderCodec()
                    .decode(BufferedLineReader.fromString(sequenceDictionaryString), file.getFile().getAbsolutePath())
                    .getSequenceDictionary();
        });
        intervals = new Lazy<>(() -> HDF5Utils.readIntervals(file, INTERVALS_GROUP_NAME));
        counts = new Lazy<>(() -> new Array2DRowRealMatrix(file.readDoubleMatrix(COUNTS_PATH)));
    }

    SampleLocatableMetadata getMetadata() {
        return new SimpleSampleLocatableMetadata(sampleName.get(), sequenceDictionary.get());
    }

    List<SimpleInterval> getIntervals() {
        return intervals.get();
    }

    /**
     * @return single-row matrix containing the counts
     */
    RealMatrix getCounts() {
        return counts.get();
    }

    /**
     * @param intervals note that no particular sort order is assumed or checked for here,
     *                  but this package-protected method should only be called by {@link SimpleCountCollection#writeHDF5},
     *                  which enforces the order specified by {@link AbstractSampleLocatableCollection}
     */
    static void write(final File outFile,
                      final SampleLocatableMetadata metadata,
                      final List<SimpleInterval> intervals,
                      final double[] counts) {
        Utils.nonNull(outFile);
        Utils.nonNull(metadata);
        Utils.nonEmpty(intervals);
        Utils.nonNull(counts);

        Utils.validateArg(intervals.size() == counts.length, "Number of intervals and counts must match.");
        Utils.validateArg(intervals.stream().distinct().count() == intervals.size(), "Intervals must all be unique.");
        Utils.validateArg(Arrays.stream(counts).noneMatch(c -> c < 0), "Counts must all be non-negative integers.");

        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            final HDF5SimpleCountCollection hdf5CountCollection = new HDF5SimpleCountCollection(file);
            hdf5CountCollection.writeSampleName(SAMPLE_NAME_PATH, metadata.getSampleName());
            hdf5CountCollection.writeSequenceDictionary(SEQUENCE_DICTIONARY_PATH, metadata.getSequenceDictionary());
            hdf5CountCollection.writeIntervals(intervals);
            hdf5CountCollection.writeCounts(counts);
        }
    }

    private <T extends SimpleInterval> void writeIntervals(final List<T> intervals) {
        HDF5Utils.writeIntervals(file, INTERVALS_GROUP_NAME, intervals);
    }

    private void writeSampleName(final String path, final String sampleName) {
        file.makeStringArray(path, sampleName);
    }

    private void writeSequenceDictionary(final String path, final SAMSequenceDictionary sequenceDictionary) {
        final StringWriter stringWriter = new StringWriter();
        new SAMTextHeaderCodec().encode(stringWriter, new SAMFileHeader(sequenceDictionary));
        file.makeStringArray(path, stringWriter.toString());
    }

    private void writeCounts(final double[] counts) {
        file.makeDoubleMatrix(COUNTS_PATH, new double[][]{counts});
    }
}
