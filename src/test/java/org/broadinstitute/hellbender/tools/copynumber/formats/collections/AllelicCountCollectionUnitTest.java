package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Unit tests for {@link AllelicCountCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicCountCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir + "copynumber/formats/collections");
    private static final File ALLELIC_COUNTS_FILE = new File(TEST_SUB_DIR, "allelic-count-collection-normal.tsv");
    private static final File ALLELIC_COUNTS_MISSING_NUCLEOTIDES_FILE = new File(TEST_SUB_DIR, "allelic-count-collection-normal-missing-nucleotides.tsv");

    private static final SampleLocatableMetadata METADATA_EXPECTED = new SimpleSampleLocatableMetadata(
            "test-sample",
            new SAMSequenceDictionary(Arrays.asList(
                    new SAMSequenceRecord("1", 20000),
                    new SAMSequenceRecord("2", 20000))));

    private static final AllelicCountCollection ALLELIC_COUNTS_EXPECTED = new AllelicCountCollection(
            METADATA_EXPECTED,
            Arrays.asList(
                    new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0, Nucleotide.G, Nucleotide.A),
                    new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4, Nucleotide.G, Nucleotide.C),
                    new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6, Nucleotide.G, Nucleotide.A),
                    new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18, Nucleotide.T, Nucleotide.A),
                    new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8, Nucleotide.C, Nucleotide.G),
                    new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8, Nucleotide.T, Nucleotide.A),
                    new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4, Nucleotide.C, Nucleotide.G),
                    new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9, Nucleotide.T, Nucleotide.A),
                    new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5, Nucleotide.G, Nucleotide.A),
                    new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0, Nucleotide.G, Nucleotide.A),
                    new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3, Nucleotide.T, Nucleotide.C)));

    private static final AllelicCountCollection ALLELIC_COUNTS_MISSING_NUCLEOTIDES_EXPECTED = new AllelicCountCollection(
            METADATA_EXPECTED,
            Arrays.asList(
                    new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0),
                    new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4),
                    new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6),
                    new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18),
                    new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8),
                    new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8),
                    new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4),
                    new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9),
                    new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5),
                    new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0),
                    new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3)));

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
        ALLELIC_COUNTS_EXPECTED.write(outputFile);
        Assert.assertTrue(FileUtils.contentEquals(outputFile, ALLELIC_COUNTS_FILE));
    }
}