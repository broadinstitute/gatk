package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberFormatsUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Unit tests for {@link AbstractSampleLocatableCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AbstractSampleLocatableCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");
    private static final File SIMPLE_LOCATABLE_COLLECTION_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection.tsv");
    private static final File SIMPLE_LOCATABLE_COLLECTION_NON_DICTIONARY_ORDER_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection-non-dictionary-order.tsv");
    private static final File SIMPLE_LOCATABLE_COLLECTION_MISSING_COLUMN_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection-missing-column.tsv");

    private static final SampleLocatableMetadata METADATA_EXPECTED = new SimpleSampleLocatableMetadata(
            "test-sample",
            new SAMSequenceDictionary(Arrays.asList(
                    new SAMSequenceRecord("1", 20000),
                    new SAMSequenceRecord("2", 20000),
                    new SAMSequenceRecord("10", 20000))));

    private static final SimpleSampleLocatableCollection SIMPLE_LOCATABLE_COLLECTION_EXPECTED = new SimpleSampleLocatableCollection(
            METADATA_EXPECTED,
            Arrays.asList(
                    new SimpleLocatable(new SimpleInterval("1", 1, 1), 1.),
                    new SimpleLocatable(new SimpleInterval("1", 2, 2), 2.),
                    new SimpleLocatable(new SimpleInterval("2", 1, 1), 3.),
                    new SimpleLocatable(new SimpleInterval("10", 1, 1), Double.NaN)));

    //simple example of a record class
    private static final class SimpleLocatable implements Locatable {
        private final SimpleInterval interval;
        private final double value;

        private SimpleLocatable(final SimpleInterval interval, final double value) {
            this.interval = Utils.nonNull(interval);
            this.value = value;
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }

        public SimpleInterval getInterval() {
            return interval;
        }

        public double getValue() {
            return value;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final SimpleLocatable that = (SimpleLocatable) o;

            return Double.compare(that.value, value) == 0 && interval.equals(that.interval);
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            result = interval.hashCode();
            temp = Double.doubleToLongBits(value);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }

        @Override
        public String toString() {
            return "SimpleLocatable{" +
                    "interval=" + interval +
                    ", value=" + value +
                    '}';
        }
    }

    //simple example of a collection class
    private static final class SimpleSampleLocatableCollection extends AbstractSampleLocatableCollection<SimpleLocatable> {
        enum SimpleLocatableTableColumn {
            CONTIG,
            START,
            END,
            VALUE;

            static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
        }
        
        private static final Function<DataLine, SimpleLocatable> SIMPLE_LOCATABLE_RECORD_FROM_DATA_LINE_DECODER = dataLine -> {
            final String contig = dataLine.get(SimpleLocatableTableColumn.CONTIG);
            final int start = dataLine.getInt(SimpleLocatableTableColumn.START);
            final int end = dataLine.getInt(SimpleLocatableTableColumn.END);
            final double value = dataLine.getDouble(SimpleLocatableTableColumn.VALUE);
            final SimpleInterval interval = new SimpleInterval(contig, start, end);
            return new SimpleLocatable(interval, value);
        };

        private static final BiConsumer<SimpleLocatable, DataLine> SIMPLE_LOCATABLE_RECORD_TO_DATA_LINE_ENCODER = (simpleLocatable, dataLine) ->
                dataLine.append(simpleLocatable.getInterval().getContig())
                        .append(simpleLocatable.getInterval().getStart())
                        .append(simpleLocatable.getInterval().getEnd())
                        .append(formatDouble(simpleLocatable.getValue()));

        private SimpleSampleLocatableCollection(final File inputFile) {
            super(inputFile, SimpleLocatableTableColumn.COLUMNS, SIMPLE_LOCATABLE_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_LOCATABLE_RECORD_TO_DATA_LINE_ENCODER);
        }

        private SimpleSampleLocatableCollection(final SampleLocatableMetadata metadata,
                                                final List<SimpleLocatable> simpleLocatables) {
            super(metadata, simpleLocatables, SimpleLocatableTableColumn.COLUMNS, SIMPLE_LOCATABLE_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_LOCATABLE_RECORD_TO_DATA_LINE_ENCODER);
        }
    }

    @Test
    public void testRead() {
        final SimpleSampleLocatableCollection simpleLocatableCollection = new SimpleSampleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_FILE);
        assertSimpleLocatableCollectionEqualsExpected(simpleLocatableCollection);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReadIntervalsNotInDictionaryOrder() {
        new SimpleSampleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_NON_DICTIONARY_ORDER_FILE);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadMissingColumn() {
        new SimpleSampleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_MISSING_COLUMN_FILE);
    }

    /**
     * Note that this will fail if {@link CopyNumberFormatsUtils#DOUBLE_FORMAT} is changed.
     */
    @Test
    public void testWrite() throws IOException {
        final File tempFile = createTempFile("test", ".tsv");
        SIMPLE_LOCATABLE_COLLECTION_EXPECTED.write(tempFile);
        SIMPLE_LOCATABLE_COLLECTION_EXPECTED.write(tempFile);   //test that file is overwritten
        Assert.assertTrue(FileUtils.contentEquals(tempFile, SIMPLE_LOCATABLE_COLLECTION_FILE));
    }

    @Test
    public void testConstructorFromListDictionarySortingOfIntervals() {
        final SimpleSampleLocatableCollection simpleLocatableCollectionExpectedUnsortedListArgument = new SimpleSampleLocatableCollection(
                METADATA_EXPECTED,
                Arrays.asList(
                        new SimpleLocatable(new SimpleInterval("1", 1, 1), 1.),
                        new SimpleLocatable(new SimpleInterval("1", 2, 2), 2.),
                        new SimpleLocatable(new SimpleInterval("10", 1, 1), Double.NaN),
                        new SimpleLocatable(new SimpleInterval("2", 1, 1), 3.)));
        assertSimpleLocatableCollectionEqualsExpected(simpleLocatableCollectionExpectedUnsortedListArgument);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalsWithDuplicates() {
        final List<SimpleLocatable> intervalsWithDuplicates = Arrays.asList(
                new SimpleLocatable(new SimpleInterval("1", 1, 1), 1.),
                new SimpleLocatable(new SimpleInterval("1", 1, 1), 1.),
                new SimpleLocatable(new SimpleInterval("2", 1, 1), 1.));
        new SimpleSampleLocatableCollection(METADATA_EXPECTED, intervalsWithDuplicates);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalsWithOverlaps() {
        final List<SimpleLocatable> intervalsWithOverlaps = Arrays.asList(
                new SimpleLocatable(new SimpleInterval("1", 1, 100), 1.),
                new SimpleLocatable(new SimpleInterval("1", 100, 200), 1.),
                new SimpleLocatable(new SimpleInterval("2", 1, 1), 1.));
        new SimpleSampleLocatableCollection(METADATA_EXPECTED, intervalsWithOverlaps);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalOutsideSequenceDictionary() {
        final List<SimpleLocatable> intervalOutsideSequenceDictionary = Collections.singletonList(
                new SimpleLocatable(new SimpleInterval("X", 1, 100), 1.));
        new SimpleSampleLocatableCollection(METADATA_EXPECTED, intervalOutsideSequenceDictionary);
    }

    private static void assertSimpleLocatableCollectionEqualsExpected(final SimpleSampleLocatableCollection simpleLocatableCollection) {
        Assert.assertEquals(simpleLocatableCollection, SIMPLE_LOCATABLE_COLLECTION_EXPECTED);
        Assert.assertEquals(simpleLocatableCollection.getMetadata(), SIMPLE_LOCATABLE_COLLECTION_EXPECTED.getMetadata());
        Assert.assertEquals(simpleLocatableCollection.getRecords(), SIMPLE_LOCATABLE_COLLECTION_EXPECTED.getRecords());
    }
}