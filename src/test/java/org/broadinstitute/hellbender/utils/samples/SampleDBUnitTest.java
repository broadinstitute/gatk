package org.broadinstitute.hellbender.utils.samples;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;


public class SampleDBUnitTest extends GATKBaseTest {
    private static SampleDBBuilder builder;
    // all the test sample files are located here
    private File testPED = new File(getToolTestDataDir() +  "testtrio.ped");

    private static final Set<Sample> testPEDSamples = new LinkedHashSet<>(Arrays.asList(
            new Sample("kid", "fam1", "dad", "mom", Sex.MALE, Affection.AFFECTED),
            new Sample("dad", "fam1", null, null, Sex.MALE, Affection.UNAFFECTED),
            new Sample("mom", "fam1", null, null, Sex.FEMALE, Affection.AFFECTED)));

    private static final Set<Sample> testPEDFamilyF2 = new LinkedHashSet<>(Arrays.asList(
            new Sample("s2", "fam2", "d2", "m2", Sex.FEMALE, Affection.AFFECTED),
            new Sample("d2", "fam2", null, null, Sex.MALE, Affection.UNKNOWN),
            new Sample("m2", "fam2", null, null, Sex.FEMALE, Affection.UNKNOWN)
    ));

    private static final Set<Sample> testPEDFamilyF3 = new LinkedHashSet<>(Arrays.asList(
            new Sample("s1", "fam3", "d1", "m1", Sex.FEMALE, Affection.AFFECTED),
            new Sample("d1", "fam3", null, null, Sex.MALE, Affection.UNKNOWN),
            new Sample("m1", "fam3", null, null, Sex.FEMALE, Affection.UNKNOWN)
    ));

    private static final Set<Sample> testSAMSamples = new LinkedHashSet<>(Arrays.asList(
            new Sample("kid", null, null, null, Sex.UNKNOWN,   Affection.UNKNOWN),
            new Sample("mom", null, null, null, Sex.UNKNOWN,   Affection.UNKNOWN),
            new Sample("dad", null, null, null, Sex.UNKNOWN,   Affection.UNKNOWN)));

    private static final Map<String, Set<Sample>> testGetFamilies = new LinkedHashMap<>();
    static {
        testGetFamilies.put("fam1", testPEDSamples);
        testGetFamilies.put("fam2", testPEDFamilyF2);
        testGetFamilies.put("fam3", testPEDFamilyF3);
    }

    private static final Set<Sample> testKidsWithParentsFamilies2 = new LinkedHashSet<>(Arrays.asList(
            new Sample("kid", "fam1", "dad", "mom", Sex.MALE,   Affection.AFFECTED),
            new Sample("kid3", "fam5", "dad2", "mom2", Sex.MALE,   Affection.AFFECTED),
            new Sample("kid2", "fam5", "dad2", "mom2", Sex.MALE,   Affection.AFFECTED)));

    private static final Set<String> testGetPartialFamiliesIds = new LinkedHashSet<>(Arrays.asList("kid", "s1"));
    private static final Map<String, Set<Sample>> testGetPartialFamilies = new LinkedHashMap<>();
    static {
        testGetPartialFamilies.put("fam1", new LinkedHashSet<>(Arrays.asList(new Sample("kid", "fam1", "dad", "mom", Sex.MALE, Affection.AFFECTED))));
        testGetPartialFamilies.put("fam3", new LinkedHashSet<>(Arrays.asList(new Sample("s1", "fam3", "d1", "m1", Sex.FEMALE, Affection.AFFECTED))));
    }

    private static final String testPEDString =
            String.format("%s%n%s%n%s",
                    "fam1 kid dad mom 1 2",
                    "fam1 dad 0   0   1 1",
                    "fam1 mom 0   0   2 2");

    private static final String testPEDMultipleFamilies =
            String.format("%s%n%s%n%s%n%s%n%s",
                    "fam1 kid dad mom 1 2",
                    "fam1 dad 0   0   1 1",
                    "fam1 mom 0   0   2 2",
                    "fam3 s1  d1  m1  2 2",
                    "fam2 s2  d2  m2  2 2");

    private static final String testPEDMultipleFamilies2 =
            String.format("%s%n%s%n%s%n%s%n%s%n%s%n%s%n%s%n%s",
                    "fam1 kid dad mom 1 2",
                    "fam1 dad 0   0   1 1",
                    "fam1 mom 0   0   2 2",
                    "fam4 kid4 dad4 0 1 2",
                    "fam4 dad4 0   0   1 1",
                    "fam5 kid2 dad2 mom2 1 2",
                    "fam5 kid3 dad2 mom2 1 2",
                    "fam5 dad2 0   0   1 1",
                    "fam5 mom2 0   0   2 2");

    private static final String testPEDStringInconsistentGender =
          "fam1 kid 0   0   2 2";

    private static final String testPEDStringConsistent =
            "fam1 kid dad   mom   1 2";

    private static final Set<Sample> testPEDSamplesAsSet =
            new LinkedHashSet<>(testPEDSamples);


    @BeforeMethod
    public void before() {
        builder = new SampleDBBuilder(PedigreeValidationType.STRICT);
    }

    @Test()
    public void loadPEDFile() {
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testPEDSamplesAsSet, db.getSamples());
    }

    @Test()
    public void loadPEDString() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDString));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testPEDSamplesAsSet, db.getSamples());
    }

    @Test()
    public void testAddSamplesFromSampleNames() {
        List<String> names = new ArrayList<>();
        testSAMSamples.forEach(s -> names.add(s.getID()));
        builder.addSamplesFromSampleNames(names);
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getSamples(), testSAMSamples);
    }

    @Test()
    public void loadDuplicateData() {
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(testPEDSamples, db.getSamples());
    }

    @Test(expectedExceptions = UserException.class)
    public void loadNonExistentFile() {
        builder.addSamplesFromPedigreeFiles(Arrays.asList(GATKBaseTest.getSafeNonExistentFile("non-existence-file.txt")));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getSamples(), testSAMSamples);
    }

    @Test(expectedExceptions = UserException.class)
    public void loadInconsistentData() {
        builder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringInconsistentGender));
        builder.getFinalSampleDB();
    }

    @Test
    public void loadAndMergeConsistentData() {
        // build a temporary DB and get the resulting sample to use for test result comparison
        SampleDBBuilder tempDBBuilder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        tempDBBuilder = tempDBBuilder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringConsistent));
        SampleDB tempFinalDB = tempDBBuilder.getFinalSampleDB();
        Sample baseKidSample = tempFinalDB.getSample("kid");

        // build a sample DB and then merge in the consistent test string
        builder = new SampleDBBuilder(PedigreeValidationType.STRICT);
        builder.addSamplesFromPedigreeFiles(Arrays.asList(testPED));
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDStringConsistent));
        SampleDB finalDB = builder.getFinalSampleDB();

        Assert.assertEquals(finalDB.getSamples().size(), 3);
        Assert.assertTrue(finalDB.getSample("kid").equals(baseKidSample));
    }

    @Test()
    public void getFamilyIDs() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getFamilyIDs(), new TreeSet<>(Arrays.asList("fam1", "fam2", "fam3")));
    }

    @Test()
    public void getFamily() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getFamily("fam1"), testPEDSamplesAsSet);
    }

    @Test()
    public void getFamilies(){
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getFamilies(),testGetFamilies);
        Assert.assertEquals(db.getFamilies(null), testGetFamilies);
        Assert.assertEquals(db.getFamilies(testGetPartialFamiliesIds),testGetPartialFamilies);
    }

    @Test()
    public void testGetFounderIds(){
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies2));
        SampleDB db = builder.getFinalSampleDB();
        Assert.assertEquals(db.getFounderIds(), new LinkedHashSet<String>(Arrays.asList("dad","mom","dad2","mom2","dad4")));
    }

    @Test()
    public void loadFamilyIDs() {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        Map<String, Set<Sample>> families = db.getFamilies();
        Assert.assertEquals(families.size(), 3);
        Assert.assertEquals(families.keySet(), new TreeSet<>(Arrays.asList("fam1", "fam2", "fam3")));

        for ( final String famID : families.keySet() ) {
            final Set<Sample> fam = families.get(famID);
            Assert.assertEquals(fam.size(), 3);
            for ( final Sample sample : fam ) {
                Assert.assertEquals(sample.getFamilyID(), famID);
            }
        }
    }

    @Test
    public void testGetTrios() throws Exception {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies));
        SampleDB db = builder.getFinalSampleDB();
        final Set<Trio> trios = db.getTrios();
        Assert.assertEquals(trios.size(), 3);
    }

    @Test
    public void testGetTriosMultipleKids() throws Exception {
        builder.addSamplesFromPedigreeStrings(Arrays.asList(testPEDMultipleFamilies2));
        SampleDB db = builder.getFinalSampleDB();
        final Set<Trio> trios = db.getTrios(false);
        Assert.assertEquals(trios.size(), 3);

        final Set<Trio> triosOneKid = db.getTrios(true);
        Assert.assertEquals(triosOneKid.size(), 1);

    }
}
