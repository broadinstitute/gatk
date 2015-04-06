package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;

import java.util.*;

/**
 * <p/>
 * This is the base test class for read filter test classes.  All read
 * filter test cases should extend from this
 * class; it sets ups a header mock up to test read filtering.
 *
 * Feel free to override non-final method to modify the behavior
 * (i.e. change how read group id are formatted, or complete a header).
 *
 * <p/>
 * You can statically determine the number of read-group involved
 * in the test by calling {@link #ReadFilterTest(int)} in you constructor.
 * <p/>
 *
 * Notice that the same header object is shared by all test and
 * it is initialized by Junit (calling {@link #beforeClass()}.
 *
 * @author Valentin Ruano Rubio
 * @date May 23, 2013
 */
public class ReadFilterTest extends BaseTest {

    private static final int DEFAULT_READ_GROUP_COUNT = 5;
    private static final int DEFAULT_READER_COUNT = 1;
    private static final String DEFAULT_READ_GROUP_PREFIX = "ReadGroup";
    private static final String DEFAULT_PLATFORM_UNIT_PREFIX = "Lane";
    private static final String DEFAULT_SAMPLE_NAME_PREFIX = "Sample";
    private static final String DEFAULT_PLATFORM_PREFIX = "Platform";
    private static final int DEFAULT_CHROMOSOME_COUNT = 1;
    private static final int DEFAULT_CHROMOSOME_START_INDEX = 1;
    private static final int DEFAULT_CHROMOSOME_SIZE = 1000;
    private static final String DEFAULT_SAM_FILE_FORMAT = "readfile-%3d.bam";

    private final int groupCount;

    private SAMFileHeader header;

    /**
     * Constructs a new read-filter test providing the number of read
     * groups in the file.
     *
     * @param groupCount number of read-group in the fictional SAM file,
     *                   must be equal or greater than 1.
     */
    protected ReadFilterTest(final int groupCount) {
        if (groupCount < 1) {
            throw new IllegalArgumentException(
                    "the read group count must at least be 1");
        }
        this.groupCount = groupCount;
    }

    /**
     * Returns the mock-up SAM file header for testing.
     *
     * @throws IllegalStateException if the header was not initialized
     *          invoking {@link #beforeClass()}
     * @return never <code>null</code>
     */
    protected final SAMFileHeader getHeader() {
        checkHeaderExists();
        return header;
    }

    /**
     * Construct a read filter test with the default number of groups
     *  ({@link #DEFAULT_READ_GROUP_COUNT}.
     */
    public ReadFilterTest() {
        this(DEFAULT_READ_GROUP_COUNT);
    }

    /**
     * Return the number of read groups involved in the test
     * @return <code>1</code> or greater.
     */
    protected final int getReadGroupCount() {
        return groupCount;
    }

    /**
     * Composes the Id for the read group given its index.
     *
     * This methods must return a unique distinct ID for each possible index and
     * it must be the same value each time it is invoked.
     *
     * @param index the index of the targeted read group in the range
     *              [1,{@link #getReadGroupCount()}]
     * @return never <code>null</code> and must be unique to each possible
     *         read group index.
     */
    protected String composeReadGroupId(final int index) {
        checkReadGroupIndex(index);
        return DEFAULT_READ_GROUP_PREFIX + index;
    }

    /**
     * Composes the Platform name for the read group given its index.
     *
     * This method must always return the same value give an index.
     *
     * @param index the index of the targeted read group in the range
     *              [1,{@link #getReadGroupCount()}]
     * @return never <code>null</code>.
     */
    protected String composePlatformName(final int index) {
        checkReadGroupIndex(index);
        return DEFAULT_PLATFORM_PREFIX + (((index-1)%2)+1);
    }


    /**
     * Composes the Platform unit name for the read group given its index.
     *
     * @param index the index of the targeted read group in the range
     *              [1,{@link #getReadGroupCount()}]
     * @return never <code>null</code>.
     */
    protected String composePlatformUnitName(final int index) {
        checkReadGroupIndex(index);
        return DEFAULT_PLATFORM_UNIT_PREFIX + (((index-1)%3)+1);
    }



    /**
     * Checks the correctness of a given read group index.
     *
     * A correct index is any value in the range [1,{@link #getReadGroupCount()}].
     *
     * @param index the target index.
     * @throws IllegalArgumentException if the input index is not correct.
     */
    protected final void checkReadGroupIndex(final int index) {
        checkIndex(index,groupCount,"read group");
    }


    private void checkIndex(final int index, final int max, CharSequence name) {
        if (index < 1 || index > max) {
            throw new IllegalArgumentException(
                    name + " index ("
                    + index
                    + ") is out of bounds [1," + max + "]");
        }
    }


    /**
     * Checks whether the header was initialized.
     *
     * @throws IllegalStateException if the header was not yet initialized.
     */
    protected final void checkHeaderExists() {
        if (header == null) {
            throw new IllegalArgumentException(
                    "header has not been initialized;"
                    + " beforeClass() was not invoked");
        }
    }

    /**
     * Checks whether the data source was initialized.
     *
     * @throws IllegalStateException if the data source was not yet initialized.
     */
    protected final void checkDataSourceExists() {
        if (header == null) {
            throw new IllegalArgumentException(
                    "data source has not been initialized;"
                            + " beforeClass() was not invoked");
        }
    }

    /**
     * Returns the ID for a read group given its index.
     *
     * @param index the index of the targeted read group in the range
     *              [1,{@link #getReadGroupCount()}]
     * @return never <code>null</code> and must be unique to each
     *              possible read group index.
     */
    protected final String getReadGroupId(final int index) {
        checkReadGroupIndex(index);
        return getHeader().getReadGroups().get(index - 1).getReadGroupId();
    }

    /**
     * Returns the platform name for a read group given its index.
     *
     * @param group the index of the targeted read group in the range
     *              [1,{@link #getReadGroupCount()}]
     * @return never <code>null</code>.
     */
    protected final String getPlatformName(final int group) {
        checkReadGroupIndex(group);
        return getHeader().getReadGroups().get(group - 1).getPlatform();
    }

    /**
     * Returns the platform unit for a read group given its index.
     *
     * @param group the index of the targeted read group in the range
     *              [1,{@link #getReadGroupCount()}]
     * @return never <code>null</code>.
     */
    protected final String getPlatformUnit(final int group) {
        checkReadGroupIndex(group);
        return getHeader().getReadGroups().get(group - 1).getPlatformUnit();
    }


    /**
     * Composes the mock up SAM file header.
     *
     * It must return an equivalent (equal) value each time it is invoked.
     *
     * @return never <code>null</code>.
     */
    protected SAMFileHeader composeHeader() {

        return ArtificialSAMUtils.createArtificialSamHeader(
                DEFAULT_CHROMOSOME_COUNT, DEFAULT_CHROMOSOME_START_INDEX,
                DEFAULT_CHROMOSOME_SIZE);
    }

    @BeforeClass
    public void beforeClass() {

        header = composeHeader();
        final List<String> readGroupIDs = new ArrayList<String>();
        final List<String> sampleNames = new ArrayList<String>();

        for (int i = 1; i <= getReadGroupCount(); i++) {
            final String readGroupId = composeReadGroupId(i);
            readGroupIDs.add(readGroupId);
            sampleNames.add(readGroupId);
        }

        ArtificialSAMUtils.createEnumeratedReadGroups(
                header, readGroupIDs, sampleNames);

        for (int i = 1; i <= getReadGroupCount(); i++) {
            final String readGroupId = readGroupIDs.get(i-1);
            final SAMReadGroupRecord groupRecord = header.getReadGroup(readGroupId);
            groupRecord.setAttribute("PL", composePlatformName(i));
            groupRecord.setAttribute("PU", composePlatformUnitName(i));
        }

    }

    @AfterClass
    public void afterClass() {
        header = null;
    }

    /**
     * Creates a read record.
     *
     * @param cigar the new record CIGAR.
     * @param group the new record group index that must be in the range \
     *              [1,{@link #getReadGroupCount()}]
     * @param reference the reference sequence index (0-based)
     * @param start the start position of the read alignment in the reference
     *              (1-based)
     * @return never <code>null</code>
     */
    protected SAMRecord createRead(final Cigar cigar, final int group, final int reference, final int start) {
        final SAMRecord record = ArtificialSAMUtils.createArtificialRead(cigar);
        record.setHeader(getHeader());
        record.setAlignmentStart(start);
        record.setReferenceIndex(reference);
        record.setAttribute(SAMTag.RG.toString(), getReadGroupId(group));
        return record;

    }
}
