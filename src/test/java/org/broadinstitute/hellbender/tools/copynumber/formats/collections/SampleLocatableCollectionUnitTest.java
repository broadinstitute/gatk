package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Unit tests for {@link SampleLocatableCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SampleLocatableCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");
    private static final File SIMPLE_LOCATABLE_COLLECTION_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection.tsv");
    private static final File SIMPLE_LOCATABLE_COLLECTION_NON_LEXICOGRAPHICAL_ORDER_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection-non-lexicographical-order.tsv");
    private static final File SIMPLE_LOCATABLE_COLLECTION_MISSING_COLUMN_FILE =
            new File(TEST_SUB_DIR, "locatable-collection-tsv-simple-locatable-collection-missing-column.tsv");
    private static final String SAMPLE_NAME_EXPECTED = "test";

    //simple example of a record class
    private static final class SimpleLocatable implements Locatable {
        private final SimpleInterval interval;
        private final int value;

        private SimpleLocatable(final SimpleInterval interval, final int value) {
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

            return value == that.value && interval.equals(that.interval);
        }

        @Override
        public int hashCode() {
            int result = interval.hashCode();
            result = 31 * result + value;
            return result;
        }
    }

    //simple example of a collection class
    private static final class SimpleSampleLocatableCollection extends SampleLocatableCollection<SimpleLocatable> {
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
            final int value = dataLine.getInt(SimpleLocatableTableColumn.VALUE);
            final SimpleInterval interval = new SimpleInterval(contig, start, end);
            return new SimpleLocatable(interval, value);
        };

        private static final BiConsumer<SimpleLocatable, DataLine> SIMPLE_LOCATABLE_RECORD_TO_DATA_LINE_ENCODER = (simpleLocatable, dataLine) ->
                dataLine.append(simpleLocatable.getInterval().getContig())
                        .append(simpleLocatable.getInterval().getStart())
                        .append(simpleLocatable.getInterval().getEnd())
                        .append(simpleLocatable.getValue());

        private SimpleSampleLocatableCollection(final File inputFile) {
            super(inputFile, SimpleLocatableTableColumn.COLUMNS, SIMPLE_LOCATABLE_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_LOCATABLE_RECORD_TO_DATA_LINE_ENCODER);
        }

        private SimpleSampleLocatableCollection(final SampleMetadata sampleMetadata,
                                                final List<SimpleLocatable> simpleLocatables) {
            super(sampleMetadata, simpleLocatables, SimpleLocatableTableColumn.COLUMNS, SIMPLE_LOCATABLE_RECORD_FROM_DATA_LINE_DECODER, SIMPLE_LOCATABLE_RECORD_TO_DATA_LINE_ENCODER);
        }
    }

    private static final SimpleSampleLocatableCollection SIMPLE_LOCATABLE_COLLECTION_EXPECTED = new SimpleSampleLocatableCollection(
            new SimpleSampleMetadata(SAMPLE_NAME_EXPECTED),
            Arrays.asList(
                    new SimpleLocatable(new SimpleInterval("1", 1, 1), 1),
                    new SimpleLocatable(new SimpleInterval("1", 2, 2), 2),
                    new SimpleLocatable(new SimpleInterval("10", 1, 1), 3),
                    new SimpleLocatable(new SimpleInterval("2", 1, 1), 4)));

    @Test
    public void testRead() {
        final SimpleSampleLocatableCollection simpleLocatableCollection = new SimpleSampleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_FILE);
        assertSimpleLocatableCollectionEqualsExpected(simpleLocatableCollection);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReadIntervalsNotInLexicographicalOrder() {
        new SimpleSampleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_NON_LEXICOGRAPHICAL_ORDER_FILE);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadMissingColumn() {
        new SimpleSampleLocatableCollection(SIMPLE_LOCATABLE_COLLECTION_MISSING_COLUMN_FILE);
    }

    @Test
    public void testWrite() throws IOException {
        final File tempFile = createTempFile("test", ".tsv");
        SIMPLE_LOCATABLE_COLLECTION_EXPECTED.write(tempFile);
        Assert.assertTrue(FileUtils.contentEquals(tempFile, SIMPLE_LOCATABLE_COLLECTION_FILE));
    }

    @Test
    public void testConstructorFromListLexicographicalSortingOfIntervals() {
        final SimpleSampleLocatableCollection simpleLocatableCollectionExpectedUnsortedListArgument = new SimpleSampleLocatableCollection(
                new SimpleSampleMetadata(SAMPLE_NAME_EXPECTED),
                Arrays.asList(
                        new SimpleLocatable(new SimpleInterval("1", 1, 1), 1),
                        new SimpleLocatable(new SimpleInterval("1", 2, 2), 2),
                        new SimpleLocatable(new SimpleInterval("2", 1, 1), 4),
                        new SimpleLocatable(new SimpleInterval("10", 1, 1), 3)));
        assertSimpleLocatableCollectionEqualsExpected(simpleLocatableCollectionExpectedUnsortedListArgument);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalsWithDuplicates() {
        final SampleMetadata sampleMetadata = new SimpleSampleMetadata("sampleName");
        final List<SimpleLocatable> intervalsWithDuplicates = Arrays.asList(
                new SimpleLocatable(new SimpleInterval("1", 1, 1), 1),
                new SimpleLocatable(new SimpleInterval("1", 1, 1), 2),
                new SimpleLocatable(new SimpleInterval("2", 1, 1), 1));
        new SimpleSampleLocatableCollection(sampleMetadata, intervalsWithDuplicates);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testIntervalsWithOverlaps() {
        final SampleMetadata sampleMetadata = new SimpleSampleMetadata("sampleName");
        final List<SimpleLocatable> intervalsWithOverlaps = Arrays.asList(
                new SimpleLocatable(new SimpleInterval("1", 1, 100), 1),
                new SimpleLocatable(new SimpleInterval("1", 100, 200), 2),
                new SimpleLocatable(new SimpleInterval("2", 1, 1), 1));
        new SimpleSampleLocatableCollection(sampleMetadata, intervalsWithOverlaps);
    }

    private static void assertSimpleLocatableCollectionEqualsExpected(final SimpleSampleLocatableCollection simpleLocatableCollection) {
        Assert.assertEquals(simpleLocatableCollection, SIMPLE_LOCATABLE_COLLECTION_EXPECTED);
        Assert.assertEquals(simpleLocatableCollection.getSampleName(), SIMPLE_LOCATABLE_COLLECTION_EXPECTED.getSampleName());
        Assert.assertEquals(simpleLocatableCollection.getRecords(), SIMPLE_LOCATABLE_COLLECTION_EXPECTED.getRecords());
    }
}