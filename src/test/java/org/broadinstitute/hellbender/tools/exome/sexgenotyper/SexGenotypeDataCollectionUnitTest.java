package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.cloud.dataflow.sdk.repackaged.com.google.common.collect.ImmutableMap;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link SexGenotypeDataCollection}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class SexGenotypeDataCollectionUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/sexgenotyper/";
    private static final File TEST_SEX_GENOTYPE_BASIC_FILE = new File(TEST_SUB_DIR, "sex_genotypes_broadies_basic.tsv");
    private static final File TEST_SEX_GENOTYPE_EXTENDED_FILE = new File(TEST_SUB_DIR, "sex_genotypes_broadies_extended.tsv");

    private static final String[] SAMPLE_NAMES = {"David", "Laura", "Geraldine", "Aaron", "Andrey", "Foo", "Yossi"};
    private static final String[] ALL_GENOTYPES = {"MALE_GIRAFFE", "FEMALE_GIRAFFE", "NEUTRAL_GIRAFFE"};
    private static final String[] SAMPLE_GENOTYPES = {"MALE_GIRAFFE", "FEMALE_GIRAFFE", "FEMALE_GIRAFFE", "MALE_GIRAFFE",
            "MALE_GIRAFFE", "NEUTRAL_GIRAFFE", "MALE_GIRAFFE"};
    private static final List<Map<String, Double>> LOG_LIKELIHOOD_MAPS = Arrays.asList(
            ImmutableMap.<String, Double>builder()
                    .put("MALE_GIRAFFE", 9.0).put("FEMALE_GIRAFFE", 0.0).put("NEUTRAL_GIRAFFE", 2.0).build(),
            ImmutableMap.<String, Double>builder()
                    .put("MALE_GIRAFFE", 5.0).put("FEMALE_GIRAFFE", 9.0).put("NEUTRAL_GIRAFFE", 1.0).build(),
            ImmutableMap.<String, Double>builder()
                    .put("MALE_GIRAFFE", 1.0).put("FEMALE_GIRAFFE", 7.0).put("NEUTRAL_GIRAFFE", 2.0).build(),
            ImmutableMap.<String, Double>builder()
                    .put("MALE_GIRAFFE", 8.0).put("FEMALE_GIRAFFE", 1.0).put("NEUTRAL_GIRAFFE", 2.0).build(),
            ImmutableMap.<String, Double>builder()
                    .put("MALE_GIRAFFE", 8.0).put("FEMALE_GIRAFFE", 1.0).put("NEUTRAL_GIRAFFE", 1.0).build(),
            ImmutableMap.<String, Double>builder()
                    .put("MALE_GIRAFFE", 1.0).put("FEMALE_GIRAFFE", 2.0).put("NEUTRAL_GIRAFFE", 8.0).build(),
            ImmutableMap.<String, Double>builder()
                    .put("MALE_GIRAFFE", 7.0).put("FEMALE_GIRAFFE", 1.0).put("NEUTRAL_GIRAFFE", 2.0).build());

    private static SexGenotypeDataCollection basicSexGenotypeCollection, extendedSexGenotypeCollection;

    @BeforeClass
    public static void makeCollections() {
        extendedSexGenotypeCollection = new SexGenotypeDataCollection();
        for (int index = 0; index < SAMPLE_NAMES.length; index++) {
            final double[] logLikelihoodArray = {LOG_LIKELIHOOD_MAPS.get(index).get("MALE_GIRAFFE"),
                    LOG_LIKELIHOOD_MAPS.get(index).get("FEMALE_GIRAFFE"), LOG_LIKELIHOOD_MAPS.get(index).get("NEUTRAL_GIRAFFE")};
            final SexGenotypeData dat = new SexGenotypeData(SAMPLE_NAMES[index], SAMPLE_GENOTYPES[index],
                    Arrays.asList(ALL_GENOTYPES), logLikelihoodArray);
            extendedSexGenotypeCollection.add(dat);
        }

        basicSexGenotypeCollection = new SexGenotypeDataCollection();
        for (int index = 0; index < SAMPLE_NAMES.length; index++) {
            final SexGenotypeData dat = new SexGenotypeData(SAMPLE_NAMES[index], SAMPLE_GENOTYPES[index], null, null);
            basicSexGenotypeCollection.add(dat);
        }
    }

    @Test
    public void testGetSampleSexGenotypeData() {
        IntStream.range(0, SAMPLE_NAMES.length).forEach(index -> {
            final String sampleName = SAMPLE_NAMES[index];
            final SexGenotypeData sexGenotypeData = extendedSexGenotypeCollection.getSampleSexGenotypeData(sampleName);
            Assert.assertEquals(sexGenotypeData, extendedSexGenotypeCollection.getSexGenotypeDataList().get(index));
            for (final String gen : ALL_GENOTYPES) {
                Assert.assertEquals(sexGenotypeData.getLogLikelihoodPerGenotype(gen),
                        LOG_LIKELIHOOD_MAPS.get(index).get(gen), 1e-12);
            }

        });
    }

    @Test
    public void testReaderBasic() {
        try {
            final SexGenotypeDataCollection col = new SexGenotypeDataCollection(TEST_SEX_GENOTYPE_BASIC_FILE);
            performBasicAsserts(col);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource files");
        }
    }

    @Test
    public void testReaderExtended() {
        try {
            final SexGenotypeDataCollection col = new SexGenotypeDataCollection(TEST_SEX_GENOTYPE_EXTENDED_FILE);
            performExtendedAsserts(col);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource files");
        }
    }

    @Test
    public void testWriteBasic() {
        final File outFile = createTempFile("foo", ".tsv");
        try {
            basicSexGenotypeCollection.write(new FileWriter(outFile));
            final SexGenotypeDataCollection reloaded = new SexGenotypeDataCollection(outFile);
            performBasicAsserts(reloaded);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outFile, "Could not create output sex genotyping file");
        }
    }

    @Test
    public void testWriteExtended() {
        final File outFile = createTempFile("foo", ".tsv");
        try {
            extendedSexGenotypeCollection.write(new FileWriter(outFile));
            final SexGenotypeDataCollection reloaded = new SexGenotypeDataCollection(outFile);
            performExtendedAsserts(reloaded);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outFile, "Could not create output sex genotyping file");
        }
    }

    private void performBasicAsserts(final SexGenotypeDataCollection col) {
        final List<SexGenotypeData> dat = col.getSexGenotypeDataList();
        Assert.assertEquals(dat.stream().map(SexGenotypeData::getSampleName).collect(Collectors.toList()),
                Arrays.asList(SAMPLE_NAMES));
        Assert.assertEquals(dat.stream().map(SexGenotypeData::getSexGenotype).collect(Collectors.toList()),
                Arrays.asList(SAMPLE_GENOTYPES));
    }

    private void performExtendedAsserts(final SexGenotypeDataCollection col) {
        final List<SexGenotypeData> dat = col.getSexGenotypeDataList();
        Assert.assertEquals(dat.stream().map(SexGenotypeData::getSampleName).collect(Collectors.toList()),
                Arrays.asList(SAMPLE_NAMES));
        Assert.assertEquals(dat.stream().map(SexGenotypeData::getSexGenotype).collect(Collectors.toList()),
                Arrays.asList(SAMPLE_GENOTYPES));
        for (int index = 0; index < SAMPLE_NAMES.length; index++) {
            final int finalIndex = index;
            final SexGenotypeData currentData = dat.stream().filter(entry -> entry.getSampleName().equals(SAMPLE_NAMES[finalIndex]))
                    .collect(Collectors.toList()).get(0);
            for (final String gen : ALL_GENOTYPES) {
                Assert.assertEquals(currentData.getLogLikelihoodPerGenotype(gen), LOG_LIKELIHOOD_MAPS.get(index).get(gen), 1e-12);
            }
        }
    }

}
