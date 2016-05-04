package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.omg.CORBA.SystemException;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link BayesianHetPulldownCalculator}.
 *
 * The fake pileups and the likelihoods are generated based on the constants defined below and for the
 * HETEROGENEOUS prior. At the moment, there are no unit test the simpler case of BALANCED prior.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class BayesianHetPulldownCalculatorUnitTest extends BaseTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    /*****************************************************************************************************************/

    private static final File SNP_FILE = new File(TEST_SUB_DIR, "common_SNP.interval_list");
    private static final File REF_FILE = new File(hg19MiniReference);
    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "normal.sorted.bam");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "tumor.sorted.bam");

    private static final int READ_DEPTH_THRESHOLD = 10;
    private static final int MINIMUM_MAPPING_QUALITY = 30;
    private static final int MINIMUM_BASE_QUALITY = 20;
    private static final int QUADRATURE_ORDER = 300;
    private static final double MIN_ABNORMAL_FRACTION = 0.5;
    private static final double MAX_ABNORMAL_FRACTION = 0.8;
    private static final int MAX_COPY_NUMBER = 4;
    private static final double ERROR_PROBABILITY_ADJUSTMENT_FACTOR = 1.0;

    private static SAMFileHeader normalHeader;
    private static SAMFileHeader tumorHeader;

    private static BayesianHetPulldownCalculator calculator;

    private static int numPileupEntries;
    private static List<Map<Nucleotide, List<BayesianHetPulldownCalculator.BaseQuality>>>
            fakePileupBaseQualities = new ArrayList<>();
    private static ArrayList<Double> fakePileupHetLogLikelihoodArray = new ArrayList<>();
    private static ArrayList<Double> fakePileupHomLogLikelihoodArray = new ArrayList<>();

    /*****************************************************************************************************************/

    private static final File FAKE_PILEUP_FILE = new File(TEST_SUB_DIR, "test_fake_generated_pileup.txt");

    public interface StringMapper <T> {
        T of(String input);
    }

    @BeforeClass
    public void initHeaders() throws IOException {
        try (final SamReader normalBamReader = SamReaderFactory.makeDefault().open(NORMAL_BAM_FILE);
             final SamReader tumorBamReader = SamReaderFactory.makeDefault().open(TUMOR_BAM_FILE)) {
            normalHeader = normalBamReader.getFileHeader();
            tumorHeader = tumorBamReader.getFileHeader();
        }
    }

    @BeforeClass
    public void initHetPulldownCalculator() {
        calculator = new BayesianHetPulldownCalculator(REF_FILE, IntervalList.fromFile(SNP_FILE),
                MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, READ_DEPTH_THRESHOLD,
                ValidationStringency.STRICT, ERROR_PROBABILITY_ADJUSTMENT_FACTOR,
                new HeterogeneousHeterozygousPileupPriorModel(MIN_ABNORMAL_FRACTION, MAX_ABNORMAL_FRACTION,
                        MAX_COPY_NUMBER, QUADRATURE_ORDER));
    }

    /**
     * load the fake fileup from file
     *
     * @throws IOException
     */
    @BeforeClass
    public void loadPileup() throws IOException {

        /* load fake pileup and likelihoods from disk */
        StringMapper<String> stringStripper = s -> s.substring(s.indexOf(':') + 1).replaceAll("\\s+", "");
        StringMapper<Double> parseToDouble = s -> Double.parseDouble(stringStripper.of(s));
        StringMapper<ArrayList<Double>> parseToDoubleArray = s -> {
            String[] tokenizedLine = stringStripper.of(s).split(",");
            ArrayList<Double> errorList = new ArrayList<>();
            if (tokenizedLine.length >= 1 && !tokenizedLine[0].equals("")) {
                errorList.addAll(Arrays.asList(tokenizedLine).stream()
                        .map(Double::parseDouble)
                        .collect(Collectors.toList()));
            }
            return errorList;
        };

        Scanner reader = new Scanner(new FileInputStream(FAKE_PILEUP_FILE));
        while (reader.hasNextLine()) {

            Map<Nucleotide, List<BayesianHetPulldownCalculator.BaseQuality>> baseQualities = new HashMap<>();

            for (Nucleotide base : new Nucleotide[]{Nucleotide.A, Nucleotide.C, Nucleotide.T, Nucleotide.G}) {

                /* load base read error list */
                ArrayList<Double> readErrorList = parseToDoubleArray.of(reader.nextLine());

                /* set the mapping error probability to 1e-6 */
                ArrayList<Double> mappingErrorList = new ArrayList<>();
                mappingErrorList.addAll(IntStream.range(0, readErrorList.size())
                        .mapToDouble(i -> 1e-6).boxed().collect(Collectors.toList()));

                /* contruct the BaseQuality list */
                List<BayesianHetPulldownCalculator.BaseQuality> baseQualityList = new ArrayList<>();
                baseQualityList.addAll(IntStream.range(0, readErrorList.size()).mapToObj(i -> new BayesianHetPulldownCalculator.BaseQuality(
                        readErrorList.get(i), mappingErrorList.get(i))).collect(Collectors.toList()));

                baseQualities.put(base, baseQualityList);
            }

            fakePileupBaseQualities.add(baseQualities);
            fakePileupHetLogLikelihoodArray.add(parseToDouble.of(reader.nextLine()));
            fakePileupHomLogLikelihoodArray.add(parseToDouble.of(reader.nextLine()));
        }

        numPileupEntries = fakePileupBaseQualities.size();
    }

    @DataProvider(name = "inputTestGetHomLogLikelihood")
    public Object[][] inputTestGetHomLogLikelihood() {

        Nucleotide alleleRef = Nucleotide.A;
        Nucleotide alleleAlt = Nucleotide.T;
        double homRefPrior = 0.5;

        Object[][] out = new Object[numPileupEntries][];
        for (int i=0; i<numPileupEntries; i++) {
            out[i] = new Object[]{
                    fakePileupBaseQualities.get(i),
                    alleleRef,
                    alleleAlt,
                    homRefPrior,
                    fakePileupHomLogLikelihoodArray.get(i)
            };
        }

        return out;
    }

    @Test(dataProvider = "inputTestGetHomLogLikelihood")
    public void testGetHomLogLikelihood(final Map<Nucleotide, List<BayesianHetPulldownCalculator.BaseQuality>> baseQualities,
                                        final Nucleotide alleleRef, final Nucleotide alleleAlt,
                                        final double homRefPrior, final double expectedHomLogLikelihood) {
        for (int i=0; i<numPileupEntries; i++) {
            Assert.assertEquals(calculator.getHomLogLikelihood(baseQualities, alleleRef, alleleAlt,
                    homRefPrior), expectedHomLogLikelihood, 1e-4);
        }
    }

    @DataProvider(name = "inputTestGetHetLogLikelihood")
    public Object[][] inputTestGetHetLogLikelihood() {

        Nucleotide alleleRef = Nucleotide.A;
        Nucleotide alleleAlt = Nucleotide.T;

        Object[][] out = new Object[numPileupEntries][];
        for (int i=0; i<numPileupEntries; i++) {
            out[i] = new Object[]{
                    fakePileupBaseQualities.get(i),
                    alleleRef,
                    alleleAlt,
                    fakePileupHetLogLikelihoodArray.get(i)
            };
        }

        return out;
    }

    @Test(dataProvider = "inputTestGetHetLogLikelihood")
    public void testGetHetLogLikelihood(final Map<Nucleotide, List<BayesianHetPulldownCalculator.BaseQuality>> baseQualities,
                                        final Nucleotide alleleRef, final Nucleotide alleleAlt,
                                        final double expectedHetLogLikelihood) {
        for (int i=0; i<numPileupEntries; i++) {
            Assert.assertEquals(calculator.getHetLogLikelihood(baseQualities, alleleRef, alleleAlt),
                    expectedHetLogLikelihood, 1e-2);
        }
    }

    @Test
    public void testInferAlleleAltFromPileup() {

        Map<Nucleotide, List<BayesianHetPulldownCalculator.BaseQuality>> baseQualities = new HashMap<>();
        BayesianHetPulldownCalculator.BaseQuality genericBaseQuality = new BayesianHetPulldownCalculator.BaseQuality(60.0, 60.0);

        List<BayesianHetPulldownCalculator.BaseQuality> genericBaseQualityList_0 = new ArrayList<>();
        List<BayesianHetPulldownCalculator.BaseQuality> genericBaseQualityList_1 = new ArrayList<>();
        List<BayesianHetPulldownCalculator.BaseQuality> genericBaseQualityList_2 = new ArrayList<>();
        List<BayesianHetPulldownCalculator.BaseQuality> genericBaseQualityList_3 = new ArrayList<>();
        List<BayesianHetPulldownCalculator.BaseQuality> genericBaseQualityList_4 = new ArrayList<>();

        /* only the base counts matter; so we can use a generic BaseQuality for all lists */

        genericBaseQualityList_1.add(genericBaseQuality);

        genericBaseQualityList_2.add(genericBaseQuality);
        genericBaseQualityList_2.add(genericBaseQuality);

        genericBaseQualityList_3.add(genericBaseQuality);
        genericBaseQualityList_3.add(genericBaseQuality);
        genericBaseQualityList_3.add(genericBaseQuality);

        genericBaseQualityList_4.add(genericBaseQuality);
        genericBaseQualityList_4.add(genericBaseQuality);
        genericBaseQualityList_4.add(genericBaseQuality);
        genericBaseQualityList_4.add(genericBaseQuality);

        /* test 1 */
        baseQualities.put(Nucleotide.A, genericBaseQualityList_4);
        baseQualities.put(Nucleotide.C, genericBaseQualityList_2);
        baseQualities.put(Nucleotide.T, genericBaseQualityList_1);
        baseQualities.put(Nucleotide.G, genericBaseQualityList_1);

        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.A), Nucleotide.C);
        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.C), Nucleotide.A);
        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.T), Nucleotide.A);
        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.G), Nucleotide.A);

        /* test 2 */
        baseQualities.put(Nucleotide.A, genericBaseQualityList_0);
        baseQualities.put(Nucleotide.C, genericBaseQualityList_3);
        baseQualities.put(Nucleotide.T, genericBaseQualityList_0);
        baseQualities.put(Nucleotide.G, genericBaseQualityList_0);

        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.A), Nucleotide.C);
        /* in this case, A is chosen simply because sorting places A right after C; in theory, T and G are equally
           valid "guesses" */
        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.C), Nucleotide.A);
        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.T), Nucleotide.C);
        Assert.assertEquals(BayesianHetPulldownCalculator.inferAltFromPileup(baseQualities, Nucleotide.G), Nucleotide.C);

    }

    @Test
    public void testGetHetPulldown() {

        Pulldown testPulldown, expectedPulldown;
        File tempFile;

        /* write Pulldown to file while checking the results */
        try {

            /* test 1: normal, loose threshold */
            testPulldown = calculator.getHetPulldown(NORMAL_BAM_FILE, 2);
            tempFile = File.createTempFile("testPulldownNormalLoose", ".txt");
            testPulldown.write(tempFile, AllelicCountTableColumns.AllelicCountTableVerbosity.FULL);

            expectedPulldown = new Pulldown(normalHeader);
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4,
                    Nucleotide.G, Nucleotide.A, 11, 18.38));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                    Nucleotide.G, Nucleotide.T, 14, 28.84));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                    Nucleotide.T, Nucleotide.G, 17, 39.39));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                    Nucleotide.T, Nucleotide.G, 15, 28.23));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5,
                    Nucleotide.G, Nucleotide.C, 11, 24.54));

            Assert.assertEquals(new Pulldown(tempFile, normalHeader), expectedPulldown);

            /* test 2: normal, tight threshold */
            testPulldown = calculator.getHetPulldown(NORMAL_BAM_FILE, 12);
            tempFile = File.createTempFile("testPulldownNormalTight", ".txt");
            testPulldown.write(tempFile, AllelicCountTableColumns.AllelicCountTableVerbosity.FULL);

            expectedPulldown = new Pulldown(normalHeader);
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                    Nucleotide.G, Nucleotide.T, 14, 28.84));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                    Nucleotide.T, Nucleotide.G, 17, 39.39));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                    Nucleotide.T, Nucleotide.G, 15, 28.23));

            Assert.assertEquals(new Pulldown(tempFile, normalHeader), expectedPulldown);

            /* test 3: tumor, loose threshold */
            testPulldown = calculator.getHetPulldown(TUMOR_BAM_FILE, 2);
            tempFile = File.createTempFile("testPulldownTumorLoose", ".txt");
            testPulldown.write(tempFile, AllelicCountTableColumns.AllelicCountTableVerbosity.FULL);

            expectedPulldown = new Pulldown(tumorHeader);
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 11522, 11522), 7, 4,
                    Nucleotide.G, Nucleotide.A, 11, 15.63));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 12098, 12098), 8, 6,
                    Nucleotide.G, Nucleotide.T, 14, 24.72));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                    Nucleotide.T, Nucleotide.G, 17, 33.89));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("2", 14689, 14689), 6, 9,
                    Nucleotide.T, Nucleotide.G, 15, 24.11));
            expectedPulldown.add(new AllelicCount(new SimpleInterval("2", 14982, 14982), 6, 5,
                    Nucleotide.G, Nucleotide.C, 11, 21.10));

            Assert.assertEquals(new Pulldown(tempFile, normalHeader), expectedPulldown);

            /* test 4: tumor, tight threshold */
            testPulldown = calculator.getHetPulldown(TUMOR_BAM_FILE, 12);
            tempFile = File.createTempFile("testPulldownTumorTight", ".txt");
            testPulldown.write(tempFile, AllelicCountTableColumns.AllelicCountTableVerbosity.FULL);

            expectedPulldown = new Pulldown(tumorHeader);
            expectedPulldown.add(new AllelicCount(new SimpleInterval("1", 14630, 14630), 9, 8,
                    Nucleotide.T, Nucleotide.G, 17, 33.89));

            Assert.assertEquals(new Pulldown(tempFile, normalHeader), expectedPulldown);

        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not write pulldown to to file.", e);
        }

    }
}
