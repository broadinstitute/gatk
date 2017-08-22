package org.broadinstitute.hellbender.utils.samples;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.StringReader;
import java.util.*;


/**
 * UnitTest for PedReader
 *
 * @author Mark DePristo
 * @since 2011
 */
public class PedReaderUnitTest extends GATKBaseTest {
    private static Logger logger = LogManager.getLogger(PedReaderUnitTest.class);

    private class PedReaderTest extends TestDataProvider {
        public String fileContents;
        public List<Sample> expectedSamples;
        EnumSet<PedReader.MissingPedField> missing;

        private PedReaderTest(final String name, final List<Sample> expectedSamples, final String fileContents) {
            super(PedReaderTest.class, name);
            this.fileContents = fileContents;
            this.expectedSamples = expectedSamples;
        }
    }

//     Family ID
//     Individual ID
//     Paternal ID
//     Maternal ID
//     Sex (1=male; 2=female; other=unknown)
//     Phenotype
//
//     -9 missing
//     0 missing
//     1 unaffected
//     2 affected

    @DataProvider(name = "readerTest")
    public Object[][] createPEDFiles() {
        new PedReaderTest("singleRecordMale",
                Arrays.asList(new Sample("kid", "fam1", null, null, Sex.MALE, Affection.UNAFFECTED)),
                "fam1 kid 0 0 1 1");

        new PedReaderTest("singleRecordFemale",
                Arrays.asList(new Sample("kid", "fam1", null, null, Sex.FEMALE, Affection.UNAFFECTED)),
                "fam1 kid 0 0 2 1");

        new PedReaderTest("singleRecordMissingGender",
                Arrays.asList(new Sample("kid", "fam1", null, null, Sex.UNKNOWN, Affection.UNKNOWN)),
                "fam1 kid 0 0 0 0");

        // Affection
        new PedReaderTest("singleRecordAffected",
                Arrays.asList(new Sample("kid", "fam1", null, null, Sex.MALE, Affection.AFFECTED)),
                "fam1 kid 0 0 1 2");

        new PedReaderTest("singleRecordUnaffected",
                Arrays.asList(new Sample("kid", "fam1", null, null, Sex.MALE, Affection.UNAFFECTED)),
                "fam1 kid 0 0 1 1");

        new PedReaderTest("singleRecordMissingAffection-9",
                Arrays.asList(new Sample("kid", "fam1", null, null, Sex.MALE, Affection.UNKNOWN)),
                "fam1 kid 0 0 1 -9");

        new PedReaderTest("singleRecordMissingAffection0",
                Arrays.asList(new Sample("kid", "fam1", null, null, Sex.MALE, Affection.UNKNOWN)),
                "fam1 kid 0 0 1 0");

        new PedReaderTest("multipleUnrelated",
                Arrays.asList(
                        new Sample("s1", "fam1", null, null, Sex.MALE,   Affection.UNAFFECTED),
                        new Sample("s2", "fam2", null, null, Sex.FEMALE, Affection.AFFECTED)),
                String.format("%s%n%s",
                        "fam1 s1 0 0 1 1",
                        "fam2 s2 0 0 2 2"));

        new PedReaderTest("multipleUnrelatedExtraLine",
                Arrays.asList(
                        new Sample("s1", "fam1", null, null, Sex.MALE,   Affection.UNAFFECTED),
                        new Sample("s2", "fam2", null, null, Sex.FEMALE, Affection.AFFECTED)),
                String.format("%s%n%s%n  %n", // note extra newlines and whitespace
                        "fam1 s1 0 0 1 1",
                        "fam2 s2 0 0 2 2"));

        new PedReaderTest("explicitTrio",
                Arrays.asList(
                        new Sample("kid", "fam1", "dad", "mom", Sex.MALE,   Affection.AFFECTED),
                        new Sample("dad", "fam1", null, null,   Sex.MALE,   Affection.UNAFFECTED),
                        new Sample("mom", "fam1", null, null,   Sex.FEMALE, Affection.AFFECTED)),
                String.format("%s%n%s%n%s",
                        "fam1 kid dad mom 1 2",
                        "fam1 dad 0   0   1 1",
                        "fam1 mom 0   0   2 2"));

        new PedReaderTest("implicitTrio",
                Arrays.asList(
                        new Sample("kid", "fam1", "dad", "mom", Sex.MALE,   Affection.AFFECTED),
                        new Sample("dad", "fam1", null, null,   Sex.MALE,   Affection.UNKNOWN),
                        new Sample("mom", "fam1", null, null,   Sex.FEMALE, Affection.UNKNOWN)),
                "fam1 kid dad mom 1 2");

        new PedReaderTest("partialTrio",
                Arrays.asList(
                        new Sample("kid", "fam1", "dad", "mom", Sex.MALE,   Affection.AFFECTED),
                        new Sample("dad", "fam1", null, null,   Sex.MALE,   Affection.UNAFFECTED),
                        new Sample("mom", "fam1", null, null,   Sex.FEMALE, Affection.UNKNOWN)),
                String.format("%s%n%s",
                        "fam1 kid dad mom 1 2",
                        "fam1 dad 0   0   1 1"));

        new PedReaderTest("bigPedigree",
                Arrays.asList(
                        new Sample("kid", "fam1", "dad",       "mom",      Sex.MALE,   Affection.AFFECTED),
                        new Sample("dad", "fam1", "granddad1", "grandma1", Sex.MALE,   Affection.UNAFFECTED),
                        new Sample("granddad1", "fam1", null, null,        Sex.MALE,   Affection.UNKNOWN),
                        new Sample("grandma1",  "fam1", null, null,        Sex.FEMALE,   Affection.UNKNOWN),
                        new Sample("mom", "fam1", "granddad2", "grandma2", Sex.FEMALE, Affection.AFFECTED),
                        new Sample("granddad2", "fam1", null, null,        Sex.MALE,   Affection.UNKNOWN),
                        new Sample("grandma2",  "fam1", null, null,        Sex.FEMALE,   Affection.UNKNOWN)),
                String.format("%s%n%s%n%s",
                        "fam1 kid dad       mom      1 2",
                        "fam1 dad granddad1 grandma1 1 1",
                        "fam1 mom granddad2 grandma2 2 2"));

        return PedReaderTest.getTests(PedReaderTest.class);
    }

    private static final void runTest(PedReaderTest test, String myFileContents, EnumSet<PedReader.MissingPedField> missing) {
        logger.warn("Test " + test);
        PedReader reader = new PedReader();
        SampleDB sampleDB = new SampleDB();
        List<Sample> readSamples = reader.parse(myFileContents, missing, sampleDB);
        Assert.assertEquals(new LinkedHashSet<Sample>(test.expectedSamples), new LinkedHashSet<Sample>(readSamples));
    }

    @Test(enabled = true, dataProvider = "readerTest")
    public void testPedReader(PedReaderTest test) {
        runTest(test, test.fileContents, EnumSet.noneOf(PedReader.MissingPedField.class));
    }

    @Test(enabled = true, dataProvider = "readerTest")
    public void testPedReaderWithComments(PedReaderTest test) {
        runTest(test, String.format("#comment%n%s", test.fileContents), EnumSet.noneOf(PedReader.MissingPedField.class));
    }

    @Test(enabled = true, dataProvider = "readerTest")
    public void testPedReaderWithSemicolons(PedReaderTest test) {
        runTest(test,
                test.fileContents.replace(String.format("%n"), ";"),
                EnumSet.noneOf(PedReader.MissingPedField.class));
    }

    // -----------------------------------------------------------------
    // missing format field tests
    // -----------------------------------------------------------------

    private class PedReaderTestMissing extends TestDataProvider {
        public EnumSet<PedReader.MissingPedField> missingDesc;
        public EnumSet<PedReader.Field> missingFields;
        public final String fileContents;
        public Sample expected;


        private PedReaderTestMissing(final String name, final String fileContents,
                                     EnumSet<PedReader.MissingPedField> missingDesc,
                                     EnumSet<PedReader.Field> missingFields,
                                     final Sample expected) {
            super(PedReaderTestMissing.class, name);
            this.fileContents = fileContents;
            this.missingDesc = missingDesc;
            this.missingFields = missingFields;
            this.expected = expected;
        }
    }

    @DataProvider(name = "readerTestMissing")
    public Object[][] createPEDFilesWithMissing() {
        new PedReaderTestMissing("missingFam",
                "fam1 kid dad mom 1 2",
                EnumSet.of(PedReader.MissingPedField.NO_FAMILY_ID),
                EnumSet.of(PedReader.Field.FAMILY_ID),
                new Sample("kid", null, "dad", "mom", Sex.MALE, Affection.AFFECTED));

        new PedReaderTestMissing("missingParents",
                "fam1 kid dad mom 1 2",
                EnumSet.of(PedReader.MissingPedField.NO_PARENTS),
                EnumSet.of(PedReader.Field.PATERNAL_ID, PedReader.Field.MATERNAL_ID),
                new Sample("kid", "fam1", null, null, Sex.MALE, Affection.AFFECTED));

        new PedReaderTestMissing("missingSex",
                "fam1 kid dad mom 1 2",
                EnumSet.of(PedReader.MissingPedField.NO_SEX),
                EnumSet.of(PedReader.Field.GENDER),
                new Sample("kid", "fam1", "dad", "mom", Sex.UNKNOWN, Affection.AFFECTED));

        new PedReaderTestMissing("missingPhenotype",
                "fam1 kid dad mom 1 2",
                EnumSet.of(PedReader.MissingPedField.NO_PHENOTYPE),
                EnumSet.of(PedReader.Field.PHENOTYPE),
                new Sample("kid", "fam1", "dad", "mom", Sex.MALE, Affection.UNKNOWN));

        new PedReaderTestMissing("missingEverythingButGender",
                "fam1 kid dad mom 1 2",
                EnumSet.of(PedReader.MissingPedField.NO_PHENOTYPE, PedReader.MissingPedField.NO_PARENTS, PedReader.MissingPedField.NO_FAMILY_ID),
                EnumSet.of(PedReader.Field.FAMILY_ID, PedReader.Field.PATERNAL_ID, PedReader.Field.MATERNAL_ID, PedReader.Field.PHENOTYPE),
                new Sample("kid", null, null, null, Sex.MALE, Affection.UNKNOWN));


        return PedReaderTestMissing.getTests(PedReaderTestMissing.class);
    }

    @Test(enabled = true, dataProvider = "readerTestMissing")
    public void testPedReaderWithMissing(PedReaderTestMissing test) {
        final String contents = sliceContents(test.missingFields, test.fileContents);
        logger.warn("Test " + test);
        PedReader reader = new PedReader();
        SampleDB sampleDB = new SampleDB();
        reader.parse(new StringReader(contents), test.missingDesc, sampleDB);
        final Sample missingSample = sampleDB.getSample("kid");
        Assert.assertEquals(test.expected, missingSample, "Missing field value not expected value for " + test);
    }

    private final static String sliceContents(EnumSet<PedReader.Field> missingFieldsSet, String full) {
        List<String> parts = new ArrayList<String>(Arrays.asList(full.split("\\s+")));
        final List<PedReader.Field> missingFields = new ArrayList<PedReader.Field>(missingFieldsSet);
        Collections.reverse(missingFields);
        for ( PedReader.Field field : missingFields )
            parts.remove(field.ordinal());
        return Utils.join("\t", parts);
    }

    // -----------------------------------------------------------------
    // parsing tags
    // -----------------------------------------------------------------

    private class PedReaderTestTagParsing extends TestDataProvider {
        public EnumSet<PedReader.MissingPedField> expected;
        public final List<String> tags;

        private PedReaderTestTagParsing(final List<String> tags, EnumSet<PedReader.MissingPedField> missingDesc) {
            super(PedReaderTestTagParsing.class);
            this.tags = tags;
            this.expected = missingDesc;
        }
    }

    @DataProvider(name = "readerTestTagParsing")
    public Object[][] createReaderTestTagParsing() {
        new PedReaderTestTagParsing(
                Collections.<String>emptyList(),
                EnumSet.noneOf(PedReader.MissingPedField.class));

        new PedReaderTestTagParsing(
                Arrays.asList("NO_FAMILY_ID"),
                EnumSet.of(PedReader.MissingPedField.NO_FAMILY_ID));

        new PedReaderTestTagParsing(
                Arrays.asList("NO_PARENTS"),
                EnumSet.of(PedReader.MissingPedField.NO_PARENTS));

        new PedReaderTestTagParsing(
                Arrays.asList("NO_PHENOTYPE"),
                EnumSet.of(PedReader.MissingPedField.NO_PHENOTYPE));

        new PedReaderTestTagParsing(
                Arrays.asList("NO_SEX"),
                EnumSet.of(PedReader.MissingPedField.NO_SEX));

        new PedReaderTestTagParsing(
                Arrays.asList("NO_SEX", "NO_PHENOTYPE"),
                EnumSet.of(PedReader.MissingPedField.NO_SEX, PedReader.MissingPedField.NO_PHENOTYPE));

        new PedReaderTestTagParsing(
                Arrays.asList("NO_SEX", "NO_PHENOTYPE", "NO_PARENTS"),
                EnumSet.of(PedReader.MissingPedField.NO_SEX, PedReader.MissingPedField.NO_PHENOTYPE, PedReader.MissingPedField.NO_PARENTS));

        return PedReaderTestTagParsing.getTests(PedReaderTestTagParsing.class);
    }

    @Test(enabled = true, dataProvider = "readerTestTagParsing")
    public void testPedReaderTagParsing(PedReaderTestTagParsing test) {
        EnumSet<PedReader.MissingPedField> parsed = PedReader.parseMissingFieldTags("test", test.tags);
        Assert.assertEquals(test.expected, parsed, "Failed to properly parse tags " + test.tags);
    }

    @Test(enabled = true, expectedExceptions = CommandLineException.class)
    public void testPedReaderTagParsing1() {
        EnumSet<PedReader.MissingPedField> parsed = PedReader.parseMissingFieldTags("test", Arrays.asList("XXX"));
    }

    @Test(enabled = true, expectedExceptions = CommandLineException.class)
    public void testPedReaderTagParsing2() {
        EnumSet<PedReader.MissingPedField> parsed = PedReader.parseMissingFieldTags("test", Arrays.asList("NO_SEX", "XXX"));
    }
}