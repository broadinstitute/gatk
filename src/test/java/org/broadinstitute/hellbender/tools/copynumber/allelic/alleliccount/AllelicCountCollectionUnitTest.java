package org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Unit tests for {@link AllelicCountCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCountCollectionUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/allelic";
    private static final File ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "allelic-count-collection-normal.tsv");
    private static final File ALLELIC_COUNTS_MISSING_NUCLEOTIDES_FILE = new File(TEST_SUB_DIR, "allelic-count-collection-normal-missing-nucleotides.tsv");
    private static final AllelicCountCollection ALLELIC_COUNTS_EXPECTED = new AllelicCountCollection();
    private static final AllelicCountCollection ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED = new AllelicCountCollection();

    static {
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0, Nucleotide.G, Nucleotide.A));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4, Nucleotide.G, Nucleotide.C));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6, Nucleotide.G, Nucleotide.A));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18, Nucleotide.T, Nucleotide.A));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8, Nucleotide.C, Nucleotide.G));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8, Nucleotide.T, Nucleotide.A));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4, Nucleotide.C, Nucleotide.G));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9, Nucleotide.T, Nucleotide.A));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5, Nucleotide.G, Nucleotide.A));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0, Nucleotide.G, Nucleotide.A));
        ALLELIC_COUNTS_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3, Nucleotide.T, Nucleotide.C));

        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0));
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3));
    }

    @Test
    public void testRead() {
        final AllelicCountCollection allelicCounts = new AllelicCountCollection(ALLELIC_COUNTS_FILE);
        Assert.assertEquals(allelicCounts, ALLELIC_COUNTS_EXPECTED);                        //nucleotides used for equality check if available
        Assert.assertEquals(allelicCounts, ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED);    //nucleotides not used for equality check if not available
    }

    //files must have all columns (including nucleotides) in order to be read in as an AllelicCountCollection
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadMissingNucleotides() {
        new AllelicCountCollection(ALLELIC_COUNTS_MISSING_NUCLEOTIDES_FILE);
    }

    @Test
    public void testWrite() throws IOException {
        final File outputFile = createTempFile("allelic-count-collection-test-output", ".tsv");
        final AllelicCountCollection allelicCounts = new AllelicCountCollection(ALLELIC_COUNTS_FILE);
        allelicCounts.write(outputFile);
        Assert.assertTrue(FileUtils.contentEquals(outputFile, ALLELIC_COUNTS_FILE));
    }

    //AllelicCountCollections must have all fields (including nucleotides) in order to be written to file
    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testWriteMissingFields() {
        final File outputFile = createTempFile("allelic-count-collection-test-output", ".tsv");
        ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED.write(outputFile);
    }
}