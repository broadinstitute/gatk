package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Integration test for {@link CollectAllelicCounts}.  Uses BAM and SNP files generated from hg19mini using wgsim.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CollectAllelicCountsIntegrationTest extends CommandLineProgramTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/copynumber/allelic";
    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-normal.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-tumor.bam");
    private static final File SITES_FILE = new File(TEST_SUB_DIR, "collect-allelic-counts-sites.interval_list");
    private static final File REF_FILE = new File(hg19MiniReference);

    @DataProvider(name = "testData")
    public Object[][] testData() throws IOException {
        //counts from IGV with minMQ = 30 and minBQ = 20
        final AllelicCountCollection normalCountsExpected = new AllelicCountCollection();
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0, Nucleotide.G, Nucleotide.N));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4, Nucleotide.G, Nucleotide.A));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6, Nucleotide.G, Nucleotide.T));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18, Nucleotide.T, Nucleotide.C));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8, Nucleotide.C, Nucleotide.A));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8, Nucleotide.T, Nucleotide.G));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4, Nucleotide.C, Nucleotide.A));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9, Nucleotide.T, Nucleotide.G));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5, Nucleotide.G, Nucleotide.C));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0, Nucleotide.G, Nucleotide.N));
        normalCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3, Nucleotide.T, Nucleotide.A));

        final AllelicCountCollection tumorCountsExpected = new AllelicCountCollection();
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0, Nucleotide.G, Nucleotide.N));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4, Nucleotide.G, Nucleotide.A));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6, Nucleotide.G, Nucleotide.T));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 17, Nucleotide.T, Nucleotide.C));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8, Nucleotide.C, Nucleotide.A));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8, Nucleotide.T, Nucleotide.G));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 3, Nucleotide.C, Nucleotide.A));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9, Nucleotide.T, Nucleotide.G));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5, Nucleotide.G, Nucleotide.C));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0, Nucleotide.G, Nucleotide.N));
        tumorCountsExpected.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3, Nucleotide.T, Nucleotide.A));

        //counts from IGV with minMQ = 30 and minBQ = 20, without nucleotides
        final AllelicCountCollection normalCountsExpectedWithoutNucleotides = new AllelicCountCollection();
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 18));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 4));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0));
        normalCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3));

        final AllelicCountCollection tumorCountsExpectedWithoutNucleotides = new AllelicCountCollection();
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 10736, 10736), 0, 0));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 12444, 12444), 0, 17));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 13059, 13059), 0, 8));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("1", 15204, 15204), 4, 3));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 15110, 15110), 6, 0));
        tumorCountsExpectedWithoutNucleotides.add(new AllelicCount(new SimpleInterval("2", 15629, 15629), 5, 3));

        return new Object[][]{
                {NORMAL_BAM_FILE, normalCountsExpected},
                {TUMOR_BAM_FILE, tumorCountsExpected},
                {NORMAL_BAM_FILE, normalCountsExpectedWithoutNucleotides},
                {TUMOR_BAM_FILE, tumorCountsExpectedWithoutNucleotides}
        };
    }

    @Test(dataProvider = "testData")
    public void test(final File inputBAMFile,
                     final AllelicCountCollection countsExpected) {
        final File outputFile = createTempFile("collect-allelic-counts-test-output", ".tsv");
        final String[] arguments = {
                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, inputBAMFile.getAbsolutePath(),
                "-L", SITES_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.REFERENCE_SHORT_NAME, REF_FILE.getAbsolutePath(),
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, outputFile.getAbsolutePath()
        };
        runCommandLine(arguments);
        final AllelicCountCollection countsResult = new AllelicCountCollection(outputFile);
        Assert.assertEquals(countsExpected, countsResult);
    }
}