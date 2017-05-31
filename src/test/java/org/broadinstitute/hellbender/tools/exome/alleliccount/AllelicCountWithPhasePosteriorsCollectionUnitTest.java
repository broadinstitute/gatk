package org.broadinstitute.hellbender.tools.exome.alleliccount;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Unit tests for {@link AllelicCountWithPhasePosteriorsCollection}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountWithPhasePosteriorsCollectionUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    private static final File BASIC_SNPS_FILE = new File(TEST_SUB_DIR, "snps-basic-with-phase-posteriors.tsv");
    private static final File INTERMEDIATE_SNPS_FILE = new File(TEST_SUB_DIR, "snps-intermediate-with-phase-posteriors.tsv");
    private static final File FULL_SNPS_FILE = new File(TEST_SUB_DIR, "snps-full-with-phase-posteriors.tsv");
    private static final File SNPS_WITH_MISSING_COLUMN_FILE = new File(TEST_SUB_DIR, "snps-with-missing-column.tsv");
    
    private static final double logProbA = Math.log(0.4);
    private static final double logProbB = Math.log(0.5);
    private static final double logProbC = Math.log(0.1);
    private static final double logProbD = Math.log(Double.MIN_VALUE);

    @Test(expectedExceptions = UserException.class)
    public void testReadWithMissingColumns() {
        new AllelicCountWithPhasePosteriorsCollection(SNPS_WITH_MISSING_COLUMN_FILE);
    }

    @Test
    public void testReadAndWriteBasicCollection() {
        final File tempFile = createTempFile("allelic-count-with-phase-posteriors-collection-test", "tsv");
        final AllelicCountWithPhasePosteriorsCollection counts = new AllelicCountWithPhasePosteriorsCollection(BASIC_SNPS_FILE);
        counts.write(tempFile, AllelicCountTableColumn.AllelicCountTableVerbosity.BASIC);

        final AllelicCountWithPhasePosteriorsCollection countsWritten = new AllelicCountWithPhasePosteriorsCollection(tempFile);

        final AllelicCountWithPhasePosteriorsCollection countsExpected = new AllelicCountWithPhasePosteriorsCollection();
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 881918, 881918), 14, 21), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 909238, 909238), 13, 11), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 934940, 934940), 14, 0), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 949608, 949608), 20, 14), logProbD, logProbB, logProbB));

        Assert.assertEquals(countsWritten, countsExpected);
    }

    @Test
    public void testReadAndWriteIntermediateCollection() {
        final File tempFile = createTempFile("allelic-count-with-phase-posteriors-collection-test", "tsv");
        final AllelicCountWithPhasePosteriorsCollection counts = new AllelicCountWithPhasePosteriorsCollection(INTERMEDIATE_SNPS_FILE);
        counts.write(tempFile, AllelicCountTableColumn.AllelicCountTableVerbosity.INTERMEDIATE);

        final AllelicCountWithPhasePosteriorsCollection countsWritten = new AllelicCountWithPhasePosteriorsCollection(tempFile);

        final AllelicCountWithPhasePosteriorsCollection countsExpected = new AllelicCountWithPhasePosteriorsCollection();
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 881918, 881918), 14, 21, Nucleotide.A, Nucleotide.C, 35), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 909238, 909238), 13, 11, Nucleotide.T, Nucleotide.G, 26), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 934940, 934940), 14, 0, Nucleotide.A, Nucleotide.G, 17), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 949608, 949608), 20, 14, Nucleotide.G, Nucleotide.T, 34), logProbD, logProbB, logProbB));

        Assert.assertEquals(countsWritten, countsExpected);
    }

    @Test
    public void testReadAndWriteFullCollection() {
        final File tempFile = createTempFile("allelic-count-with-phase-posteriors-collection-test", "tsv");
        final AllelicCountWithPhasePosteriorsCollection counts = new AllelicCountWithPhasePosteriorsCollection(FULL_SNPS_FILE);
        counts.write(tempFile, AllelicCountTableColumn.AllelicCountTableVerbosity.FULL);

        final AllelicCountWithPhasePosteriorsCollection countsWritten = new AllelicCountWithPhasePosteriorsCollection(tempFile);

        final AllelicCountWithPhasePosteriorsCollection countsExpected = new AllelicCountWithPhasePosteriorsCollection();
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 881918, 881918), 14, 21, Nucleotide.A, Nucleotide.C, 35, 12.0), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 909238, 909238), 13, 11, Nucleotide.T, Nucleotide.G, 26, 120.5), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 934940, 934940), 14, 0, Nucleotide.A, Nucleotide.G, 17, -16.1), logProbA, logProbB, logProbC));
        countsExpected.add(new AllelicCountWithPhasePosteriors(new AllelicCount(new SimpleInterval("1", 949608, 949608), 20, 14, Nucleotide.G, Nucleotide.T, 34, 54.3), logProbD, logProbB, logProbB));

        Assert.assertEquals(countsWritten, countsExpected);
    }
}
