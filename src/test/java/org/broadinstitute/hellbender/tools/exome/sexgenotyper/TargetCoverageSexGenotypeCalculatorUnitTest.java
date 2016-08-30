package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollection;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link TargetCoverageSexGenotypeCalculator}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class TargetCoverageSexGenotypeCalculatorUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/sexgenotyper/";
    private static final File TEST_CONTIG_PLOIDY_ANNOTS_FILE = new File(TEST_SUB_DIR, "contig_annots.tsv");
    private static final File TEST_RCC_FILE = new File(TEST_SUB_DIR, "sex_genotyper_rcc_trunc.tsv");
    private static final File TEST_SEX_GENOTYPE_FILE = new File(TEST_SUB_DIR, "sex_genotypes_agilent_trunc.tsv");

    private static final double DEFAULT_MAPPING_ERROR_PROBABILITY = 1e-6;

    private static TargetCoverageSexGenotypeCalculator genotyper;

    private static final String[] GENOTYPES = new String[] {"SEX_XX", "SEX_XY"};
    private static final int AUTOSOME_PLOIDY = 2;
    private static final int SEX_XX_PLOIDY_ON_X = 2;
    private static final int SEX_XX_PLOIDY_ON_Y = 0;
    private static final int SEX_XY_PLOIDY_ON_X = 1;
    private static final int SEX_XY_PLOIDY_ON_Y = 1;

    @BeforeSuite @Override
    public void setTestVerbosity() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
    }

    @BeforeClass
    public static void initSexGenotyper() {
        final List<ContigGermlinePloidyAnnotation> contigPloidyAnnotsList;
        final ReadCountCollection readCounts;
        try {
            contigPloidyAnnotsList =
                    ContigGermlinePloidyAnnotationTableReader.readContigGermlinePloidyAnnotationsFromFile(TEST_CONTIG_PLOIDY_ANNOTS_FILE);
            readCounts = ReadCountCollectionUtils.parse(TEST_RCC_FILE);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource files");
        }
        genotyper = new TargetCoverageSexGenotypeCalculator(readCounts, contigPloidyAnnotsList, DEFAULT_MAPPING_ERROR_PROBABILITY);
    }

    /**
     * Note: this test relies on the content of the test resource files
     */
    @Test
    public void testGenotyperInitialization() {
        final Set<String> ploidyTags = genotyper.getSexGenotypeIdentifiers();

        final int[] autosomalTargetPloidies = genotyper.getAutosomalTargetGermlinePloidies();
        final Map<String, int[]> allosomalTargetPloidies = genotyper.getAllosomalTargetGermlinePloidiesMap();

        /* assert all ploidy tags are present */
        Assert.assertTrue(ploidyTags.equals(Arrays.stream(GENOTYPES).collect(Collectors.toSet())));

        /* assert ploidy on autosomal targets */
        ArrayAsserts.assertArrayEquals(
                IntStream.range(0, autosomalTargetPloidies.length).map(i -> AUTOSOME_PLOIDY).toArray(),
                autosomalTargetPloidies);

        /* assert ploidy on allosomal targets */
        int[] expected;
        final List<Target> allosomalTargetList = genotyper.getAllosomalTargetList();
        expected = allosomalTargetList.stream().mapToInt(t -> t.getContig().equals("X") ? SEX_XX_PLOIDY_ON_X :
                SEX_XX_PLOIDY_ON_Y).toArray();
        ArrayAsserts.assertArrayEquals(expected, allosomalTargetPloidies.get("SEX_XX"));
        expected = allosomalTargetList.stream().mapToInt(t -> t.getContig().equals("X") ? SEX_XY_PLOIDY_ON_X :
                SEX_XY_PLOIDY_ON_Y).toArray();
        ArrayAsserts.assertArrayEquals(expected, allosomalTargetPloidies.get("SEX_XY"));

    }

    /**
     * Note: this test relies on the content of the test resource files
     */
    @Test
    public void testInferSexGenotypes() {
        final SexGenotypeDataCollection result = genotyper.inferSexGenotypes();
        try {
            final SexGenotypeDataCollection expected = new SexGenotypeDataCollection(TEST_SEX_GENOTYPE_FILE);
            Assert.assertEquals(result, expected);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource files");
        }
    }

}
