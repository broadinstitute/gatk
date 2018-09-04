package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
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
    //note to developers: repeat the column headers in Javadoc so that they are viewable when linked
    /**
     * CONTIG, START, END, COUNT
     */
    enum SimpleCountTableColumn {
        CONTIG,
        START,
        END,
        COUNT;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
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

    public static SimpleCountCollection read(final File file) {
        IOUtils.canReadFile(file);
        if (IOUtils.isHDF5File(file.toPath())) {
            return readHDF5(new HDF5File(file));
        } else {
            return readTSV(file);
        }
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

    public void writeHDF5(final File file) {
        Utils.nonNull(file);
        HDF5SimpleCountCollection.write(file, getMetadata(), getIntervals(), getCounts());
    }

    public double[] getCounts() {
        return getRecords().stream().mapToDouble(SimpleCount::getCount).toArray();
    }
}
