package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.ImmutableList;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.copynumber.SimpleCountCodec;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple data structure to pass and read/write a List of {@link SimpleCount} objects.
 * Supports both TSV and HDF5.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SimpleCountCollection extends AbstractSampleLocatableCollection<SimpleCount> {
    private static final int DEFAULT_FEATURE_QUERY_LOOKAHEAD_IN_BP = 1_000_000;

    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, COUNT
     *
     * Note: Unlike the package-private enums in other collection classes, this enum and its
     * {@link TableColumnCollection} are public so that they can be accessed by {@link SimpleCountCodec},
     * which must be in org.broadinstitute.hellbender.utils.codecs to be discovered as a codec.
     */
    public enum SimpleCountTableColumn {
        CONTIG,
        START,
        END,
        COUNT;

        public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }
    
    private static final Function<DataLine, SimpleCount> SIMPLE_COUNT_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
        final String contig = dataLine.get(SimpleCountTableColumn.CONTIG);
        final int start = dataLine.getInt(SimpleCountTableColumn.START);
        final int end = dataLine.getInt(SimpleCountTableColumn.END);
        final int count = dataLine.getInt(SimpleCountTableColumn.COUNT);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new SimpleCount(interval, count);
    };

    private static final BiConsumer<SimpleCount, DataLine> SIMPLE_COUNT_RECORD_TO_DATA_LINE_ENCODER = (simpleCount, dataLine) ->
            dataLine.append(simpleCount.getInterval().getContig())
                    .append(simpleCount.getInterval().getStart())
                    .append(simpleCount.getInterval().getEnd())
                    .append(simpleCount.getCount());

    private SimpleCountCollection(final File inputFile) {
        super(inputFile, SimpleCountCollection.SimpleCountTableColumn.COLUMNS, SIMPLE_COUNT_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_COUNT_RECORD_TO_DATA_LINE_ENCODER);
    }

    public SimpleCountCollection(final SampleLocatableMetadata metadata,
                                 final List<SimpleCount> simpleCounts) {
        super(metadata, simpleCounts, SimpleCountCollection.SimpleCountTableColumn.COLUMNS, SIMPLE_COUNT_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_COUNT_RECORD_TO_DATA_LINE_ENCODER);
    }

    /**
     * Read all counts from a file (HDF5 or TSV).
     */
    public static SimpleCountCollection read(final File file) {
        IOUtils.canReadFile(file);
        return readAndSubset(file, null);  //specifying a null intervalSubset returns all counts
    }

    /**
     * From a file (HDF5 or TSV), subset only the counts with intervals coinciding with intervals from a given list.
     * The list may contain intervals that do not coincide with any count intervals.
     * Unlike {@link #readOverlappingSubsetFromGCS(String, List)}, this method first reads and constructs a {@link SimpleCountCollection}
     * using the entire file, and then creates and returns a second {@link SimpleCountCollection} containing only the
     * requested subset.  This method also does not subset count intervals that overlap but do not strictly coincide
     * with intervals in the given list, and does not require that intervals in the given list do not overlap each other.
     * @param intervalSubset    if {@code null} or empty, all counts will be returned
     */
    public static SimpleCountCollection readAndSubset(final File file,
                                                      final Set<SimpleInterval> intervalSubset) {
        IOUtils.canReadFile(file);
        final SimpleCountCollection simpleCounts = IOUtils.isHDF5File(file.toPath())
                ? readHDF5(new HDF5File(file))
                : readTSV(file);
        if (intervalSubset == null || intervalSubset.isEmpty()) {
            return simpleCounts;
        }
        return new SimpleCountCollection(
                simpleCounts.getMetadata(),
                simpleCounts.getRecords().stream()
                        .filter(c -> intervalSubset.contains(c.getInterval()))
                        .collect(Collectors.toList()));
    }

    private static SimpleCountCollection readTSV(final File file) {
        IOUtils.canReadFile(file);
        return new SimpleCountCollection(file);
    }

    private static SimpleCountCollection readHDF5(final HDF5File file) {
        Utils.nonNull(file);
        final HDF5SimpleCountCollection hdf5CountCollection = new HDF5SimpleCountCollection(file);
        final SampleLocatableMetadata metadata = hdf5CountCollection.getMetadata();
        final List<SimpleInterval> intervals = hdf5CountCollection.getIntervals();
        final double[] counts = hdf5CountCollection.getCounts().getRow(0);
        final List<SimpleCount> simpleCounts = IntStream.range(0, intervals.size())
                .mapToObj(i -> new SimpleCount(intervals.get(i), (int) counts[i]))
                .collect(Collectors.toList());
        return new SimpleCountCollection(metadata, simpleCounts);
    }

    /**
     * Read all counts from a Google Cloud Storage URL.
     * A corresponding index for the counts file must also be present.
     */
    public static SimpleCountCollection readFromGCS(final String path) {
        IOUtils.assertFileIsReadable(IOUtils.getPath(path));
        Utils.validate(BucketUtils.isGcsUrl(path), "Read-count path must be a Google Cloud Storage URL.");
        Utils.validate(new SimpleCountCodec().canDecode(path), String.format(
                "Read-count file extension must be one of the following: [%s]",
                String.join(",", SimpleCountCodec.SIMPLE_COUNT_CODEC_EXTENSIONS)));
        return readOverlappingSubsetFromGCS(path, null);  //specifying a null intervalSubset returns all counts
    }

    /**
     * From a Google Cloud Storage URL, subset only the counts with intervals overlapping with intervals from a given list.
     * A corresponding index for the counts file must also be present.
     * The given list may contain intervals that do not overlap with any count intervals.
     * Unlike {@link #readAndSubset(File, Set)}, this method checks for overlaps, rather than requiring the intervals
     * to be strictly coincident.  Calling code can use {@link IntervalUtils#getIntervalsWithFlanks} to create a merged
     * list of the original intervals desired to be strictly coincident; this merged list can then be used with this method.
     * @param overlapIntervals    if {@code null} or empty, all counts will be returned; must be sorted and non-overlapping otherwise
     */
    public static SimpleCountCollection readOverlappingSubsetFromGCS(final String path,
                                                                     final List<SimpleInterval> overlapIntervals) {
        IOUtils.assertFileIsReadable(IOUtils.getPath(path));
        Utils.validate(BucketUtils.isGcsUrl(path), "Read-count path must be a Google Cloud Storage URL.");
        Utils.validate(new SimpleCountCodec().canDecode(path), String.format(
                "Read-count file extension must be one of the following: [%s]",
                String.join(",", SimpleCountCodec.SIMPLE_COUNT_CODEC_EXTENSIONS)));
        final FeatureDataSource<SimpleCount> simpleCountsFeatureDataSource = new FeatureDataSource<>(
                path,
                path,
                DEFAULT_FEATURE_QUERY_LOOKAHEAD_IN_BP,
                SimpleCount.class,
                ConfigFactory.getInstance().getGATKConfig().cloudPrefetchBuffer(),
                ConfigFactory.getInstance().getGATKConfig().cloudIndexPrefetchBuffer());
        final SampleLocatableMetadata metadata = (SampleLocatableMetadata) simpleCountsFeatureDataSource.getHeader();
        if (overlapIntervals != null) {
            CopyNumberArgumentValidationUtils.validateIntervals(overlapIntervals, metadata.getSequenceDictionary());
            CopyNumberArgumentValidationUtils.validateNonOverlappingIntervals(overlapIntervals, metadata.getSequenceDictionary());
        }
        simpleCountsFeatureDataSource.setIntervalsForTraversal(overlapIntervals);
        final List<SimpleCount> simpleCounts = ImmutableList.copyOf(simpleCountsFeatureDataSource.iterator());
        return new SimpleCountCollection(metadata, simpleCounts);
    }

    public void writeHDF5(final File file) {
        Utils.nonNull(file);
        HDF5SimpleCountCollection.write(file, getMetadata(), getIntervals(), getCounts());
    }

    public double[] getCounts() {
        return getRecords().stream().mapToDouble(SimpleCount::getCount).toArray();
    }
}
