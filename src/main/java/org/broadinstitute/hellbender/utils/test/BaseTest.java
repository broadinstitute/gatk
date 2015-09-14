package org.broadinstitute.hellbender.utils.test;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.BeforeClass;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;


/**
 * This is the base test class for all of our test cases.  All test cases should extend from this
 * class; it sets up the logger, and resolves the location of directories that we rely on.
 */
@SuppressWarnings("unchecked")
public abstract class BaseTest {

    static {
        // set properties for the local Spark runner
        System.setProperty("dataflow.spark.test.reuseSparkContext", "true");
        SparkContextFactory.enableTestSparkContext();
    }

    public static final Logger logger = LogManager.getLogger("org.broadinstitute.gatk");

    private static final String CURRENT_DIRECTORY = System.getProperty("user.dir");
    public static final String gatkDirectory = System.getProperty("gatkdir", CURRENT_DIRECTORY) + "/";

    private static final String publicTestDirRelative = "src/test/resources/";
    public static final String publicTestDir = new File(gatkDirectory, publicTestDirRelative).getAbsolutePath() + "/";
    public static final String publicTestDirRoot = publicTestDir.replace(publicTestDirRelative, "");

    /**
     * LARGE FILES FOR TESTING (MANAGED BY GIT LFS)
     */
    public static final String largeFileTestDir = new File(publicTestDir, "large").getAbsolutePath() + "/";

    // All of chromosomes 20 and 21 from the b37 reference
    public static final String b37_reference_20_21 = largeFileTestDir + "human_g1k_v37.20.21.fasta";

    // ~600,000 reads from chromosomes 20 and 21 of an NA12878 WGS bam aligned to b37, plus ~50,000 unmapped reads
    public static final String NA12878_20_21_WGS_bam = largeFileTestDir + "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam";

    // Variants from a DBSNP 138 VCF overlapping the reads in NA12878_20_21_WGS_bam
    public static final String dbsnp_138_b37_20_21_vcf = largeFileTestDir + "dbsnp_138.b37.20.21.vcf";

    /**
     * END OF LARGE FILES FOR TESTING
     */

    public static final String hg19_chr1_1M_Reference = publicTestDir + "Homo_sapiens_assembly19_chr1_1M.fasta";
    public static final String hg19_chr1_1M_dbSNP = publicTestDir + "Homo_sapiens_assembly19.dbsnp135.chr1_1M.exome_intervals.vcf";

    // the following file has been modified such that the first chromosome length is 1M; this is sometimes
    // required due to sequence dictionary validation, since a reference FASTA with only 1M bases is used
    public static final String hg19_chr1_1M_dbSNP_modified = publicTestDir + "HSA19.dbsnp135.chr1_1M.exome_intervals.modified.vcf";

    public static final String hg19_chr1_1M_exampleVCF = publicTestDir + "joint_calling.chr1_1M.1kg_samples.10samples.noINFO.vcf";
    public static final String hg19MiniReference = publicTestDir + "hg19mini.fasta";

    public static final String exampleFASTA = publicTestDir + "exampleFASTA.fasta";
    public static final String exampleReference = hg19MiniReference;
    public static final String hg19MiniIntervalFile = publicTestDir + "hg19mini.interval_list";

    public CachingIndexedFastaSequenceFile hg19ReferenceReader;
    public GenomeLocParser hg19GenomeLocParser;

    // used to seed the genome loc parser with a sequence dictionary
    protected SAMFileHeader hg19Header;

    /**
     *  name of the google cloud project that stores the data and will run the code, probably broad-dsde-dev
     *  @return HELLBENDER_TEST_PROJECT env. var if defined, throws otherwise.
     */
    public static String getDataflowTestProject() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_PROJECT");
    }

    /**
     * API key for HELLBENDER_TEST_PROJECT
     * @return HELLBENDER_TEST_APIKEY env. var if defined, throws otherwise.
     */
    public static String getDataflowTestApiKey() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_APIKEY");
    }

    /**
     * A writeable folder on the project's GCS, where java files will be staged for execution
     * @return HELLBENDER_TEST_STAGING env. var if defined, throws otherwise.
     */
    public static String getDataflowTestStaging() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_STAGING");
    }

    /**
     *  A GCS path where the test inputs are stored
     *  @return HELLBENDER_TEST_INPUTS env. var if defined, throws otherwise.
     */
    public static String getDataflowTestInputPath() {
        return getNonNullEnvironmentVariable("HELLBENDER_TEST_INPUTS");
    }

    private static String getNonNullEnvironmentVariable(String envVarName) {
        String value = System.getenv(envVarName);
        if (null == value) {
            throw new UserException("For this test, please define environment variable \""+envVarName+"\"");
        }
        return value;
    }

    /**
     * @return an arguments List containing the apiKey, project, and staging arguments
     *         populated from environment variables as defined in {@link #getDataflowTestApiKey},
     *         {@link #getDataflowTestProject}, and {@link #getDataflowTestStaging}, suitable
     *         for use in a hellbender command line.
     */
    public static List<String> getStandardDataflowArgumentsFromEnvironment() {
        return Arrays.asList("--apiKey", getDataflowTestApiKey(),
                             "--project", getDataflowTestProject(),
                             "--staging", getDataflowTestStaging());
    }

    /**
     * Gets a PipelineOptions object containing our API key as specified in the HELLBENDER_TEST_APIKEY
     * environment variable. Useful for tests that need to access data in a GCS bucket via the
     * methods in the {@link org.broadinstitute.hellbender.utils.dataflow.BucketUtils} class,
     * but don't need to run an actual dataflow pipeline.
     *
     * @return a PipelineOptions object containing our API key
     */
    public static GCSOptions getAuthenticatedPipelineOptions() {
        return BucketUtils.getAuthenticatedGCSOptions(getDataflowTestApiKey());
    }

    @BeforeClass
    public void initGenomeLocParser() throws FileNotFoundException {
        hg19ReferenceReader = new CachingIndexedFastaSequenceFile(new File(hg19MiniReference));
        hg19Header = new SAMFileHeader();
        hg19Header.setSequenceDictionary(hg19ReferenceReader.getSequenceDictionary());
        hg19GenomeLocParser = new GenomeLocParser(hg19ReferenceReader);
    }

    protected List<GenomeLoc> intervalStringsToGenomeLocs( String... intervals) {
        return intervalStringsToGenomeLocs(Arrays.asList(intervals));
    }

    protected List<GenomeLoc> intervalStringsToGenomeLocs( List<String> intervals ) {
        List<GenomeLoc> locs = new ArrayList<>();
        for (String interval: intervals)
            locs.add(hg19GenomeLocParser.parseGenomeLoc(interval));
        return Collections.unmodifiableList(locs);
    }

    /**
     * Returns the location of the resource directory for the command line program.
     */
    public String getToolTestDataDir(){
        return "src/test/resources/" +this.getClass().getPackage().getName().replace(".","/") +"/" + getTestedClassName() + "/";
    }

    /**
     * Returns the name of the class being tested.
     * The default implementation takes the simple name of the test class and removes the trailing "Test".
     * Override if needed.
     */
    public String getTestedClassName(){
        if (getClass().getSimpleName().contains("IntegrationTest"))
            return getClass().getSimpleName().replaceAll("IntegrationTest$", "");
        else if (getClass().getSimpleName().contains("UnitTest"))
            return getClass().getSimpleName().replaceAll("UnitTest$", "");
        else
            return getClass().getSimpleName().replaceAll("Test$", "");
    }

    /**
     * Simple generic utility class to creating TestNG data providers:
     *
     * 1: inherit this class, as in
     *
     *      private class SummarizeDifferenceTest extends TestDataProvider {
     *         public SummarizeDifferenceTest() {
     *           super(SummarizeDifferenceTest.class);
     *         }
     *         ...
     *      }
     *
     * Provide a reference to your class to the TestDataProvider constructor.
     *
     * 2: Create instances of your subclass.  Return from it the call to getTests, providing
     * the class type of your test
     *
     * <code>
     * {@literal @}DataProvider(name = "summaries")
     * public Object[][] createSummaries() {
     *   new SummarizeDifferenceTest().addDiff("A", "A").addSummary("A:2");
     *   new SummarizeDifferenceTest().addDiff("A", "B").addSummary("A:1", "B:1");
     *   return SummarizeDifferenceTest.getTests(SummarizeDifferenceTest.class);
     * }
     * </code>
     *
     * This class magically tracks created objects of this
     */
    public static class TestDataProvider {
        private static final Map<Class<?>, List<Object>> tests = new HashMap<>();
        protected String name;

        /**
         * Create a new TestDataProvider instance bound to the class variable C
         */
        public TestDataProvider(Class<?> c, String name) {
            if ( ! tests.containsKey(c) )
                tests.put(c, new ArrayList<>());
            tests.get(c).add(this);
            this.name = name;
        }

        public TestDataProvider(Class<?> c) {
            this(c, "");
        }

        public void setName(final String name) {
            this.name = name;
        }

        /**
         * Return all of the data providers in the form expected by TestNG of type class C
         */
        public static Object[][] getTests(Class<?> c) {
            List<Object[]> params2 = new ArrayList<>();
            for ( Object x : tests.get(c) ) params2.add(new Object[]{x});
            return params2.toArray(new Object[][]{});
        }

        @Override
        public String toString() {
            return "TestDataProvider("+name+")";
        }
    }

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     *
     * This will also mark the corresponding Tribble/Tabix/BAM indices matching the temp file for deletion.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(final String name, final String extension) {
        return IOUtils.createTempFile(name, extension);
    }

    /**
     * Log this message so that it shows up inline during output as well as in html reports
     */
    public static void log(final String message) {
        Reporter.log(message, true);
    }

    private static final double DEFAULT_FLOAT_TOLERANCE = 1e-1;


    /**
     * Checks whether two double array contain the same values or not.
     * @param actual actual produced array.
     * @param expected expected array.
     * @param tolerance maximum difference between double value to be consider equivalent.
     */
    protected static void assertEqualsDoubleArray(final double[] actual, final double[] expected, final double tolerance) {
        if (expected == null)
            Assert.assertNull(actual);
        else {
            Assert.assertNotNull(actual);
            Assert.assertEquals(actual.length,expected.length,"array length");
        }
        for (int i = 0; i < actual.length; i++)
            Assert.assertEquals(actual[i],expected[i],tolerance,"array position " + i);
    }

    public static void assertEqualsDoubleSmart(final Object actual, final Double expected, final double tolerance) {
        Assert.assertTrue(actual instanceof Double, "Not a double");
        assertEqualsDoubleSmart((double) (Double) actual, (double) expected, tolerance);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected) {
        assertEqualsDoubleSmart(actual, expected, DEFAULT_FLOAT_TOLERANCE);
    }

    public static <T> void assertEqualsSet(final Set<T> actual, final Set<T> expected, final String info) {
        final Set<T> actualSet = new HashSet<>(actual);
        final Set<T> expectedSet = new HashSet<>(expected);
        Assert.assertTrue(actualSet.equals(expectedSet), info); // note this is necessary due to testng bug for set comps
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance) {
        assertEqualsDoubleSmart(actual, expected, tolerance, null);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance, final String message) {
        if ( Double.isNaN(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isNaN(actual), "expected is nan, actual is not");
        else if ( Double.isInfinite(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isInfinite(actual), "expected is infinite, actual is not");
        else {
            final double delta = Math.abs(actual - expected);
            final double ratio = Math.abs(actual / expected - 1.0);
            Assert.assertTrue(delta < tolerance || ratio < tolerance, "expected = " + expected + " actual = " + actual
                    + " not within tolerance " + tolerance
                    + (message == null ? "" : "message: " + message));
        }
    }

    private static void assertAttributeEquals(final String key, final Object actual, final Object expected) {
        if ( expected instanceof Double ) {
            // must be very tolerant because doubles are being rounded to 2 sig figs
            assertEqualsDoubleSmart(actual, (Double) expected, 1e-2);
        } else
            Assert.assertEquals(actual, expected, "Attribute " + key);
    }

    public static void assertGenotypesAreEqual(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype names");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "Genotype alleles");
        Assert.assertEquals(actual.getGenotypeString(), expected.getGenotypeString(), "Genotype string");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");

        // filters are the same
        Assert.assertEquals(actual.getFilters(), expected.getFilters(), "Genotype fields");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "Genotype isFiltered");

        // inline attributes
        Assert.assertEquals(actual.getDP(), expected.getDP(), "Genotype dp");
        Assert.assertTrue(Arrays.equals(actual.getAD(), expected.getAD()));
        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "Genotype gq");
        Assert.assertEquals(actual.hasPL(), expected.hasPL(), "Genotype hasPL");
        Assert.assertEquals(actual.hasAD(), expected.hasAD(), "Genotype hasAD");
        Assert.assertEquals(actual.hasGQ(), expected.hasGQ(), "Genotype hasGQ");
        Assert.assertEquals(actual.hasDP(), expected.hasDP(), "Genotype hasDP");

        Assert.assertEquals(actual.hasLikelihoods(), expected.hasLikelihoods(), "Genotype haslikelihoods");
        Assert.assertEquals(actual.getLikelihoodsString(), expected.getLikelihoodsString(), "Genotype getlikelihoodsString");
        Assert.assertEquals(actual.getLikelihoods(), expected.getLikelihoods(), "Genotype getLikelihoods");
        Assert.assertTrue(Arrays.equals(actual.getPL(), expected.getPL()));

        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "Genotype phredScaledQual");
        assertAttributesEquals(actual.getExtendedAttributes(), expected.getExtendedAttributes());
        Assert.assertEquals(actual.isPhased(), expected.isPhased(), "Genotype isPhased");
        Assert.assertEquals(actual.getPloidy(), expected.getPloidy(), "Genotype getPloidy");
    }

    private static void assertAttributesEquals(final Map<String, Object> actual, Map<String, Object> expected) {
        final Set<String> expectedKeys = new HashSet<>(expected.keySet());

        for ( final Map.Entry<String, Object> act : actual.entrySet() ) {
            final Object actualValue = act.getValue();
            if ( expected.containsKey(act.getKey()) && expected.get(act.getKey()) != null ) {
                final Object expectedValue = expected.get(act.getKey());
                if ( expectedValue instanceof List ) {
                    final List<Object> expectedList = (List<Object>)expectedValue;
                    Assert.assertTrue(actualValue instanceof List, act.getKey() + " should be a list but isn't");
                    final List<Object> actualList = (List<Object>)actualValue;
                    Assert.assertEquals(actualList.size(), expectedList.size(), act.getKey() + " size");
                    for ( int i = 0; i < expectedList.size(); i++ )
                        assertAttributeEquals(act.getKey(), actualList.get(i), expectedList.get(i));
                } else
                    assertAttributeEquals(act.getKey(), actualValue, expectedValue);
            } else {
                // it's ok to have a binding in x -> null that's absent in y
                Assert.assertNull(actualValue, act.getKey() + " present in one but not in the other");
            }
            expectedKeys.remove(act.getKey());
        }

        // now expectedKeys contains only the keys found in expected but not in actual,
        // and they must all be null
        for ( final String missingExpected : expectedKeys ) {
            final Object value = expected.get(missingExpected);
            Assert.assertTrue(isMissing(value), "Attribute " + missingExpected + " missing in one but not in other" );
        }
    }

    private static boolean isMissing(final Object value) {
        if ( value == null ) return true;
        else if ( value.equals(VCFConstants.MISSING_VALUE_v4) ) return true;
        else if ( value instanceof List ) {
            // handles the case where all elements are null or the list is empty
            for ( final Object elt : (List)value)
                if ( elt != null )
                    return false;
            return true;
        } else
            return false;
    }


    /**
     * captures {@link System#out} while runnable is executing
     * @param runnable a code block to execute
     * @return everything written to {@link System#out} by runnable
     */
    public static String captureStdout(Runnable runnable){
        return captureSystemStream(runnable, System.out, System::setOut);
    }


    /**
     * captures {@link System#err} while runnable is executing
     * @param runnable a code block to execute
     * @return everything written to {@link System#err} by runnable
     */
    public static String captureStderr(Runnable runnable){
        return captureSystemStream(runnable, System.err, System::setErr);
    }

    private static String captureSystemStream(Runnable runnable,  PrintStream stream, Consumer<? super PrintStream> setter){
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        setter.accept(new PrintStream(out));
        try {
            runnable.run();
        } finally{
            setter.accept(stream);
        }
        return out.toString();
    }

    public static void assertContains(String actual, String expectedSubstring){
        Assert.assertTrue(actual.contains(expectedSubstring),  expectedSubstring +" was not found in " + actual+ ".");
    }

}

